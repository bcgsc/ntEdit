CPPFLAGS= -O3 -std=c++11 -fopenmp -lz 

all: ntedit

ntedit:
	g++ $(CPPFLAGS) -o ntedit ntedit.cpp 

clean: 
	rm ntedit
