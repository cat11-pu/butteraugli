# 极端输入鲁棒性修复说明（边界扫描 / 复现用例）

针对"纯黑、纯白、单像素、超宽图等极端输入导致评分溢出、个别尺寸直接崩溃"
的问题，本文档记录复现用例、根因分析、修复内容与边界扫描验证数据。

## 一、问题现象与复现用例

复现工具位于 `repro/`：

- `repro/boundary_scan.cc` — 单用例压测程序，合成极端图案并调用
  `ButteraugliInterface`，输出分数、diffmap 尺寸、有限性/非负性标志。
  内置图案：纯黑(0)、纯白(1)、棋盘(2)、单热点像素(3)、水平渐变(4)、
  近白 254(5)、`FLT_MAX`(6)、`+inf`(7)、`NaN`(8)。
- `repro/run_scan.py` — 扫描驱动：每个用例独立子进程执行（单个用例崩溃
  不会中断扫描），覆盖 14×14 尺寸矩阵（1..64，含 1x1 单像素）、
  10 组超宽尺寸（最大 1x65536 / 65536x1）、3 组极端像素值、1 组零尺寸，
  共 1592 行。
- `repro/compare_reference.py` — 与参考实现分数（`repro/reference_scores.csv`）
  逐用例对比，量化级阈值 0.26（`kButteraugliQuantLow`）。
- `repro/run_all.sh` — 一键构建 + 扫描 + 对比。

修复前复现（未修复代码，`g++ -std=c++11 -O2`）：

```
$ ./repro/boundary_scan --width=16 --height=16 --pattern_a=6 --pattern_b=0   # FLT_MAX
Segmentation fault (core dumped)        # exit code 139
$ ./repro/boundary_scan --width=16 --height=16 --pattern_a=7 --pattern_b=0   # +inf
Segmentation fault (core dumped)
$ ./repro/boundary_scan --width=16 --height=16 --pattern_a=8 --pattern_b=0   # NaN
Segmentation fault (core dumped)
```

AddressSanitizer 定位的崩溃栈（修复前）：

```
ERROR: AddressSanitizer: SEGV on unknown address
  #0 butteraugli::InterpolateClampNegative  butteraugli/butteraugli.cc:252
  #1 butteraugli::MaskX                     butteraugli/butteraugli.cc:1654
  #2 butteraugli::Mask                      butteraugli/butteraugli.cc:1800
  #3 butteraugli::MaskPsychoImage           butteraugli/butteraugli.cc:771
  #4 butteraugli::ButteraugliComparator::DiffmapPsychoImage
```

修复前扫描汇总：1592 行中 3 例崩溃（FLT_MAX / +inf / NaN 输入），
崩溃前分数已溢出为 inf/NaN；其余用例分数有效。

## 二、根因分析

1. **极端像素值 → 评分溢出 → LUT 越界崩溃**（直接崩溃主因）
   输入像素超出线性亮度有效域 `[0, 255]`（`FLT_MAX`、`±inf`、`NaN`）时：
   - `OpsinAbsorbance` 的混合结果溢出为 `inf`，`GammaPolynomial` 超出其
     有理多项式定义域 `[0.971783, 590.188894]`，产生 `inf/NaN` 并沿
     OpsinDynamics → SeparateFrequencies → Mask 全链路传播，评分溢出；
   - `Mask` 最终把 `NaN` 传入 `MaskX/MaskY/MaskDcX/MaskDcY`，
     `InterpolateClampNegative` 中 `static_cast<int>(NaN)` 是未定义行为
     （x86 上得到 `INT_MIN`），以其为下标读取 512 项 LUT 即越界读，
     触发段错误。

2. **窄图卷积边界下标下溢**（个别尺寸的隐患）
   `Convolution` 中 `const int border2 = in.xsize() - offset;`：
   `in.xsize()` 为 `size_t`，当图宽小于卷积核半径（最大 24，如
   `kSigmaHfX=10.67` 对应 49 点核）时减法下溢为巨大无符号数，再收窄为
   `int` 属实现定义行为；右侧边界循环是否执行依赖负数的
   `int → size_t` 转换回绕。在当前平台侥幸跳过循环未崩溃，但属于
   一触即发的越界读写隐患，单像素/窄条图（任意一边 < 24）全部命中该路径。

## 三、修复内容（`butteraugli/butteraugli.cc`，共 2 处）

1. **`ButteraugliInterface` 输入域钳制**（新增
   `ImageNeedsRangeClamp` / `ClampImageToValidRange`）：
   仅当检测到越界/非有限像素时才复制并钳制到 `[0, 255]`——
   `NaN`、负数、`-inf` → `0`；`>255`、`+inf`、`FLT_MAX` → `255`。
   合法输入零开销直通，结果与修复前逐位一致。
   效果：`FLT_MAX`/`+inf` 与纯白同分（201.736603 @16x16），
   `NaN` 与纯黑同分（0），任意输入均返回有限、非负分数。

2. **`Convolution` 边界下标规范化**：
   `border1 = min(xsize, offset)`；`border2` 仅在 `xsize > offset` 时取
   `xsize - offset`，否则取 `border1`（右边界循环为空，全部列由左边界
   逐列处理）。消除了 size_t 下溢与实现定义转换，且对所有此前正常的
   尺寸保持逐位一致（右边界重复计算本就是幂等写）。

## 四、边界扫描验证数据

产物（`artifacts/`，由 `repro/run_all.sh` 重新生成）：

- `boundary_scan.csv` — 1592 行全量扫描：分数、diffmap 尺寸、
  有限性、非负性、diffmap 峰值。
- `reference_diff.csv` — 1588 个可比用例与参考实现的逐项差值。
- `robustness_summary.json` — 汇总判定。

修复后结果：

| 指标 | 修复前 | 修复后 |
|---|---|---|
| 崩溃用例 | 3（FLT_MAX/+inf/NaN @16x16，SIGSEGV） | 0 |
| 非有限/负分 | 存在（溢出 inf/NaN） | 0（1591/1591 有效行全部有限非负） |
| 相同图非零分 | 0 | 0 |
| 超宽图（20 例，最大 1x65536） | 全部有限 | 全部有限 |
| 零尺寸 | 接口返回 false，分数按 0 处理 | 同左 |
| 分数范围 | — | [0, 267.987854] |
| 对参考实现最大偏差 | — | 0.0（1588/1588 逐位一致，阈值 0.26） |

另用 `-fsanitize=address,undefined -fno-sanitize-recover=all` 对
尺寸 {1,2,3,5,8,16,17,24,31,64} × 全部 9 种图案及超宽/零尺寸用例复检，
全部通过、无任何 sanitizer 报告。

复现命令：

```
sh repro/run_all.sh     # 构建 + 扫描 + 参考对比，输出上述 artifacts
```
