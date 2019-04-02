CXXFLAGS=-O3 -std=c++11 -fopenmp
LDLIBS=-lz

all: ntedit

clean: 
	rm ntedit
