CXX                = mpicxx
GCC                = g++
STD                = -std=c++20
OPTFLAGS	       = -O3
CXXFLAGS          += -Wall 
CXXFLAGS          += -fopenmp 
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

test_ff: ff
	./ff 128 8
	./ff 128 16
	./ff 1024 8
	./ff 1024 16
	./ff 4096 8
	./ff 4096 16

test_omp: omp 
	mpirun -n 4 omp 128
	mpirun -n 4 omp 512
	mpirun -n 4 omp 1024
	mpirun -n 4 omp 4096

run_ff: ff 
	sbatch scripts/script_ff.sh	

run_omp: omp
	sbatch scripts/script_omp.sh	

clean: 
	-rm ff omp