CXX                = /home/simo/j/mpi_bin/bin/mpicxx
GCC                = g++
STD                = -std=c++20
OPTFLAGS	       = -O3
CXXFLAGS          += -Wall 
INCS              += -I ~/fastflow
FF_SRC             = src_ff/
OMP_SRC            = src_omp/
BIN                = bin/

.PHONY: all clean cleanall 

all: ff omp

ff:	$(FF_SRC)main.cpp
	$(GCC) $(INCS) $(STD) $(CXXFLAGS) $(OPTFLAGS) -o $@ $< 

omp: $(OMP_SRC)main.cpp
	$(CXX) $(STD) $(CXXFLAGS) $(OPTFLAGS) -o $@ $< 

clean: 
	-rm ff omp