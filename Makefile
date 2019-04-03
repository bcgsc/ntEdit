CPPFLAGS= -O3 -std=c++11 -fopenmp -lz 

all: ntedit

ntedit:
	g++ $(CPPFLAGS) -o ntedit_link ntedit_link.cpp 

clean: 
	rm ntedit_link
