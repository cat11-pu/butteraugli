LDLIBS += -lpng -ljpeg
CXXFLAGS += -std=c++11 -I..
LINK.o = $(LINK.cc)

all: butteraugli.o butteraugli_main.o butteraugli

butteraugli: butteraugli.o butteraugli_main.o

test: region_test
	./region_test

region_test: butteraugli_region.o butteraugli_region_test.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ $(LDLIBS) -o $@

butteraugli_region.o: butteraugli.cc
	$(CXX) $(CXXFLAGS) -O2 -c butteraugli.cc -o $@

butteraugli_region_test.o: butteraugli_region_test.cc
	$(CXX) $(CXXFLAGS) -O2 -c butteraugli_region_test.cc -o $@
