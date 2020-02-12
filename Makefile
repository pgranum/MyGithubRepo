#Add this path and root's path includes
CXXFLAGS   = -I./$(INCDIR) -I$(shell root-config --incdir) 
#Set optimisation and use C++ version from root (C++11 on lxplus)
CXXFLAGS +=-O2 -g  -Wall $(shell root-config --cflags)

ROOTLIBS  := $(shell root-config --libs) -Wl,-rpath,$(ROOTSYS)/lib

LDLIBS += $$(gsl-config --libs) #links GSL library
PLOTTER = gnuplot #links gnu plot

all: MySim.exe


clean:
	-rm -f *.o MySim.exe


 MySim.exe: %.exe: 
	$(CXX) MySIm.cpp -o $@ $(MACHFLAG) $(LDFLAGS) $(CXXFLAGS) $^ $(MIDASLIBS) $(ROOTLIBS) $(LDLIBS) -lm -lz -lutil  -lpthread 
