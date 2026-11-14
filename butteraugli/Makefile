LDLIBS += -lpng -ljpeg
CXXFLAGS += -std=c++11 -I..
LINK.o = $(LINK.cc)

all: butteraugli.o butteraugli_main.o butteraugli

butteraugli: butteraugli.o butteraugli_main.o

band_test: butteraugli.o butteraugli_band_test.o
	$(LINK.cc) $^ -o $@

test: band_test
	./band_test

clean:
	rm -f butteraugli butteraugli.o butteraugli_main.o band_test \
	      butteraugli_band_test.o
