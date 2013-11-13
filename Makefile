MIR=/a
CASACORE=/a/casa/local
prefix=/a

CXXFLAGS = -Wall -g -O0 -I$(MIR)/include -I$(CASACORE)/include/casacore
LFLAGS = -L$(CASACORE)/lib -L$(MIR)/lib \
 -lcasa_casa -lcasa_tables -lcasa_measures -lcasa_ms -lcasa_scimath -lcasa_scimath_f -lmir \
 -Wl,--rpath -Wl,$(CASACORE)/lib -Wl,--rpath -Wl,$(MIR)/lib

all: mirtoms mirmsflagextract

mirtoms: mirtoms.cc Makefile
	g++ -o $@ $(CXXFLAGS) $(LFLAGS) $<

mirmsflagextract: mirmsflagextract.cc Makefile
	g++ -o $@ $(CXXFLAGS) $(LFLAGS) $<

clean:
	-rm -f mirtoms mirmsflagextract

install: mirtoms mirmsflagextract
	install -m755 $^ $(prefix)/bin
