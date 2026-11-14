cc_library(
    name = "butteraugli_lib",
    srcs = [
        "butteraugli/butteraugli.cc",
        "butteraugli/butteraugli.h",
    ],
    hdrs = [
        "butteraugli/butteraugli.h",
    ],
    copts = ["-Wno-sign-compare"],
    visibility = ["//visibility:public"],
)

cc_binary(
    name = "butteraugli",
    srcs = ["butteraugli/butteraugli_main.cc"],
    copts = ["-Wno-sign-compare"],
    visibility = ["//visibility:public"],
    deps = [
        ":butteraugli_lib",
        "@jpeg_archive//:jpeg",
        "@png_archive//:png",
    ],
)

cc_test(
    name = "butteraugli_region_test",
    srcs = ["butteraugli/butteraugli_region_test.cc"],
    copts = ["-Wno-sign-compare"],
    deps = [
        ":butteraugli_lib",
    ],
)
