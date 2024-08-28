CXX                = /home/simo/openbin/bin/mpicxx
GCC                = g++
STD                = -std=c++20
OPTFLAGS	       = -O3
CXXFLAGS          += -Wall 
INCS              += -I ~/fastflow
FF_SRC             = src_ff/
BIN                = bin/

.PHONY: all clean cleanall 

ff:	$(FF_SRC)main.cpp
	$(GCC) $(INCS) $(STD) $(CXXFLAGS) $(OPTFLAGS) -o $(BIN)$@ $< 

omp: 

run: ff
	./$(BIN)ff

clean: 
	-rm -fr bin/*