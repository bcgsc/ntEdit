CXXFLAGS=-O3 -std=c++11 -fopenmp
LDLIBS=-lz

all: ntedit

# Check the C++ source code for errors.
lint: clang-format clang-tidy

# Check the C++ source code for errors with clang-tidy.
clang-tidy:
	clang-tidy -warnings-as-errors='*' *.cpp -- -std=c++11 -x c++

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in *.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
ifneq ("$(wildcard *.fixed)","")
	exit 1
endif

check:
	cd demo
	bash runme.sh

clean: 
	rm ntedit
