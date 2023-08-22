CXXFLAGS=-O3 -std=c++11 -fopenmp
LDLIBS=-lz -lboost_iostreams

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
	./ntedit -f demo/ecoliWithMismatches001Indels0001.fa.gz -r demo/solidBF_k25.bf -d 5 -i 4 -b ntEditEcolik25Test
	diff -q ntEditEcolik25Test_changes.tsv demo/nteditk25_changes.tsv && rm ntEditEcolik25Test_changes.tsv
	diff -q ntEditEcolik25Test_edited.fa demo/nteditk25_edited.fa && rm ntEditEcolik25Test_edited.fa

clean: 
	rm ntedit
