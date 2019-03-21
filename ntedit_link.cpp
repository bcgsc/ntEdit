#define PROGRAM "ntEdit"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <getopt.h>
#include <zlib.h>
#include <unordered_map>
#include <vector>
#include <queue>
#include <string>
#include <cmath>
#include <omp.h>
#include <cassert>
#include <cerrno>
#include <unistd.h>
#include "lib/kseq.h"
#include "lib/nthash.hpp"
#include "lib/BloomFilter.hpp"

KSEQ_INIT(gzFile, gzread)

static const char VERSION_MESSAGE[] = 
	PROGRAM " Version 1.1.0\n"
	"Written by Rene Warren, Hamid Mohamadi, and Jessica Zhang.\n"
	"Copyright 2018 Canada's Michael smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] = 
	"Usage: " PROGRAM " v1.1.0\n" 
	"\n"
	"Scalable genome assembly polishing.\n"
	"\n"
	"\n"
	" Options:\n"
	"	-t,	number of threads [default=1]\n"
	"	-f,	Draft genome assembly (FASTA, Multi-FASTA, and/or gzipped compatible), REQUIRED\n"
	"	-r,	Bloom filter file (generated from ntHits), REQUIRED\n"
	"	-b,	output file prefix, OPTIONAL\n"
	"	-k,	kmer size, REQUIRED\n"
	"	-z,	minimum contig length [default=100]\n"
	"	-i,	maximum number of insertion bases to try, range 0-5, [default=4]\n"
	"	-d,	maximum number of deletions bases to try, range 0-5, [default=5]\n"
	"	-x,	k/x ratio for the number of kmers that should be missing, [default=5.000]\n"
	"	-y, 	k/y ratio for the number of editted kmers that should be present, [default=9.000]\n"
	"	-m,	mode of editing, range 0-2, [default=0]\n"
	"			0: best substitution, or first good indel\n"
	"			1: best substitution, or best indel\n"
	"			2: best edit overall (suggestion that you reduce i and d for performance)\n"
	"	-v,	verbose mode (-v 1 = yes, default = 0, no)\n"
	"\n"
	"	--help,		display this message and exit \n"
	"	--version,	output version information and exit\n"
	"Report bugs to rwarren@bcgsc.ca\n"; 

using namespace std; 

namespace opt {
	unsigned nthreads=1; 
	string draft_filename; 
	string bloom_filename;
	unsigned k; 
	unsigned h=0;
	unsigned min_contig_len=100; 
	unsigned max_insertions=4; 
	unsigned max_deletions=5;
	float edit_threshold=9.0000; 
	float missing_threshold=5.0000;
	string outfile_prefix; 
	int mode=0; 
	int verbose=0; 
}

static const char shortopts[] = "t:f:s:k:z:b:r:v:d:i:x:y:m:"; 

enum { OPT_HELP = 1, OPT_VERSION }; 

static const struct option longopts[] = {
	{"threads",	required_argument, NULL, 't'},
	{"draft_file",	required_argument, NULL, 'f'},
	{"k",	required_argument, NULL, 'k'},
	{"minimum_contig_length",	required_argument, NULL, 'z'},
	{"maximum_insertions",	required_argument, NULL, 'i'},
	{"maximum_deletions", required_argument, NULL, 'd'},
	{"edit_threshold",	required_argument, NULL, 'y'},
	{"missing_threshold",	required_argument, NULL, 'x'},
	{"bloom_filename", required_argument, NULL, 'r'},
	{"outfile_prefix", required_argument, NULL, 'b'},
	{"mode", required_argument, NULL, 'm'},
	{"verbose", required_argument, NULL, 'v'},
	{"help", no_argument, NULL, OPT_HELP},
	{"version", no_argument, NULL, OPT_VERSION},
	{NULL, 0, NULL, 0}
};

// Setting up the number of tries when for each number of base insertion
vector<int> num_tries = {0,1,5,21,85,341};

// Setting up base array
unordered_map<unsigned char, vector<unsigned char>> bases_array = {{'A', {'T', 'C', 'G'}},
							    {'T', {'A', 'C', 'G'}},
							    {'C', {'A', 'T', 'G'}},
							    {'G', {'A', 'T', 'C'}},
							    {'N', {'A', 'T', 'C', 'G'}}};

// Setting all the indel combos
unordered_map<unsigned char, vector<string>> multi_possible_bases = {{'A', { "A", "AA", "AC", "AG", "AT", "AAA", "AAC", "AAG", "AAT", "ACA", "ACC", "ACG", "ACT", "AGA", "AGC", "AGG", "AGT", "ATA", "ATC", "ATG", "ATT", "AAAA", "AAAC", "AAAG", "AAAT", "AACA", "AACC", "AACG", "AACT", "AAGA", "AAGC", "AAGG", "AAGT", "AATA", "AATC", "AATG", "AATT", "ACAA", "ACAC", "ACAG", "ACAT", "ACCA", "ACCC", "ACCG", "ACCT", "ACGA", "ACGC", "ACGG", "ACGT", "ACTA", "ACTC", "ACTG", "ACTT", "AGAA", "AGAC", "AGAG", "AGAT", "AGCA", "AGCC", "AGCG", "AGCT", "AGGA", "AGGC", "AGGG", "AGGT", "AGTA", "AGTC", "AGTG", "AGTT", "ATAA", "ATAC", "ATAG", "ATAT", "ATCA", "ATCC", "ATCG", "ATCT", "ATGA", "ATGC", "ATGG", "ATGT", "ATTA", "ATTC", "ATTG", "ATTT", "AAAAA", "AAAAC", "AAAAG", "AAAAT", "AAACA", "AAACC", "AAACG", "AAACT", "AAAGA", "AAAGC", "AAAGG", "AAAGT", "AAATA", "AAATC", "AAATG", "AAATT", "AACAA", "AACAC", "AACAG", "AACAT", "AACCA", "AACCC", "AACCG", "AACCT", "AACGA", "AACGC", "AACGG", "AACGT", "AACTA", "AACTC", "AACTG", "AACTT", "AAGAA", "AAGAC", "AAGAG", "AAGAT", "AAGCA", "AAGCC", "AAGCG", "AAGCT", "AAGGA", "AAGGC", "AAGGG", "AAGGT", "AAGTA", "AAGTC", "AAGTG", "AAGTT", "AATAA", "AATAC", "AATAG", "AATAT", "AATCA", "AATCC", "AATCG", "AATCT", "AATGA", "AATGC", "AATGG", "AATGT", "AATTA", "AATTC", "AATTG", "AATTT", "ACAAA", "ACAAC", "ACAAG", "ACAAT", "ACACA", "ACACC", "ACACG", "ACACT", "ACAGA", "ACAGC", "ACAGG", "ACAGT", "ACATA", "ACATC", "ACATG", "ACATT", "ACCAA", "ACCAC", "ACCAG", "ACCAT", "ACCCA", "ACCCC", "ACCCG", "ACCCT", "ACCGA", "ACCGC", "ACCGG", "ACCGT", "ACCTA", "ACCTC", "ACCTG", "ACCTT", "ACGAA", "ACGAC", "ACGAG", "ACGAT", "ACGCA", "ACGCC", "ACGCG", "ACGCT", "ACGGA", "ACGGC", "ACGGG", "ACGGT", "ACGTA", "ACGTC", "ACGTG", "ACGTT", "ACTAA", "ACTAC", "ACTAG", "ACTAT", "ACTCA", "ACTCC", "ACTCG", "ACTCT", "ACTGA", "ACTGC", "ACTGG", "ACTGT", "ACTTA", "ACTTC", "ACTTG", "ACTTT", "AGAAA", "AGAAC", "AGAAG", "AGAAT", "AGACA", "AGACC", "AGACG", "AGACT", "AGAGA", "AGAGC", "AGAGG", "AGAGT", "AGATA", "AGATC", "AGATG", "AGATT", "AGCAA", "AGCAC", "AGCAG", "AGCAT", "AGCCA", "AGCCC", "AGCCG", "AGCCT", "AGCGA", "AGCGC", "AGCGG", "AGCGT", "AGCTA", "AGCTC", "AGCTG", "AGCTT", "AGGAA", "AGGAC", "AGGAG", "AGGAT", "AGGCA", "AGGCC", "AGGCG", "AGGCT", "AGGGA", "AGGGC", "AGGGG", "AGGGT", "AGGTA", "AGGTC", "AGGTG", "AGGTT", "AGTAA", "AGTAC", "AGTAG", "AGTAT", "AGTCA", "AGTCC", "AGTCG", "AGTCT", "AGTGA", "AGTGC", "AGTGG", "AGTGT", "AGTTA", "AGTTC", "AGTTG", "AGTTT", "ATAAA", "ATAAC", "ATAAG", "ATAAT", "ATACA", "ATACC", "ATACG", "ATACT", "ATAGA", "ATAGC", "ATAGG", "ATAGT", "ATATA", "ATATC", "ATATG", "ATATT", "ATCAA", "ATCAC", "ATCAG", "ATCAT", "ATCCA", "ATCCC", "ATCCG", "ATCCT", "ATCGA", "ATCGC", "ATCGG", "ATCGT", "ATCTA", "ATCTC", "ATCTG", "ATCTT", "ATGAA", "ATGAC", "ATGAG", "ATGAT", "ATGCA", "ATGCC", "ATGCG", "ATGCT", "ATGGA", "ATGGC", "ATGGG", "ATGGT", "ATGTA", "ATGTC", "ATGTG", "ATGTT", "ATTAA", "ATTAC", "ATTAG", "ATTAT", "ATTCA", "ATTCC", "ATTCG", "ATTCT", "ATTGA", "ATTGC", "ATTGG", "ATTGT", "ATTTA", "ATTTC", "ATTTG", "ATTTT"}},
								{'C', { "C", "CA", "CC", "CG", "CT", "CAA", "CAC", "CAG", "CAT", "CCA", "CCC", "CCG", "CCT", "CGA", "CGC", "CGG", "CGT", "CTA", "CTC", "CTG", "CTT", "CAAA", "CAAC", "CAAG", "CAAT", "CACA", "CACC", "CACG", "CACT", "CAGA", "CAGC", "CAGG", "CAGT", "CATA", "CATC", "CATG", "CATT", "CCAA", "CCAC", "CCAG", "CCAT", "CCCA", "CCCC", "CCCG", "CCCT", "CCGA", "CCGC", "CCGG", "CCGT", "CCTA", "CCTC", "CCTG", "CCTT", "CGAA", "CGAC", "CGAG", "CGAT", "CGCA", "CGCC", "CGCG", "CGCT", "CGGA", "CGGC", "CGGG", "CGGT", "CGTA", "CGTC", "CGTG", "CGTT", "CTAA", "CTAC", "CTAG", "CTAT", "CTCA", "CTCC", "CTCG", "CTCT", "CTGA", "CTGC", "CTGG", "CTGT", "CTTA", "CTTC", "CTTG", "CTTT", "CAAAA", "CAAAC", "CAAAG", "CAAAT", "CAACA", "CAACC", "CAACG", "CAACT", "CAAGA", "CAAGC", "CAAGG", "CAAGT", "CAATA", "CAATC", "CAATG", "CAATT", "CACAA", "CACAC", "CACAG", "CACAT", "CACCA", "CACCC", "CACCG", "CACCT", "CACGA", "CACGC", "CACGG", "CACGT", "CACTA", "CACTC", "CACTG", "CACTT", "CAGAA", "CAGAC", "CAGAG", "CAGAT", "CAGCA", "CAGCC", "CAGCG", "CAGCT", "CAGGA", "CAGGC", "CAGGG", "CAGGT", "CAGTA", "CAGTC", "CAGTG", "CAGTT", "CATAA", "CATAC", "CATAG", "CATAT", "CATCA", "CATCC", "CATCG", "CATCT", "CATGA", "CATGC", "CATGG", "CATGT", "CATTA", "CATTC", "CATTG", "CATTT", "CCAAA", "CCAAC", "CCAAG", "CCAAT", "CCACA", "CCACC", "CCACG", "CCACT", "CCAGA", "CCAGC", "CCAGG", "CCAGT", "CCATA", "CCATC", "CCATG", "CCATT", "CCCAA", "CCCAC", "CCCAG", "CCCAT", "CCCCA", "CCCCC", "CCCCG", "CCCCT", "CCCGA", "CCCGC", "CCCGG", "CCCGT", "CCCTA", "CCCTC", "CCCTG", "CCCTT", "CCGAA", "CCGAC", "CCGAG", "CCGAT", "CCGCA", "CCGCC", "CCGCG", "CCGCT", "CCGGA", "CCGGC", "CCGGG", "CCGGT", "CCGTA", "CCGTC", "CCGTG", "CCGTT", "CCTAA", "CCTAC", "CCTAG", "CCTAT", "CCTCA", "CCTCC", "CCTCG", "CCTCT", "CCTGA", "CCTGC", "CCTGG", "CCTGT", "CCTTA", "CCTTC", "CCTTG", "CCTTT", "CGAAA", "CGAAC", "CGAAG", "CGAAT", "CGACA", "CGACC", "CGACG", "CGACT", "CGAGA", "CGAGC", "CGAGG", "CGAGT", "CGATA", "CGATC", "CGATG", "CGATT", "CGCAA", "CGCAC", "CGCAG", "CGCAT", "CGCCA", "CGCCC", "CGCCG", "CGCCT", "CGCGA", "CGCGC", "CGCGG", "CGCGT", "CGCTA", "CGCTC", "CGCTG", "CGCTT", "CGGAA", "CGGAC", "CGGAG", "CGGAT", "CGGCA", "CGGCC", "CGGCG", "CGGCT", "CGGGA", "CGGGC", "CGGGG", "CGGGT", "CGGTA", "CGGTC", "CGGTG", "CGGTT", "CGTAA", "CGTAC", "CGTAG", "CGTAT", "CGTCA", "CGTCC", "CGTCG", "CGTCT", "CGTGA", "CGTGC", "CGTGG", "CGTGT", "CGTTA", "CGTTC", "CGTTG", "CGTTT", "CTAAA", "CTAAC", "CTAAG", "CTAAT", "CTACA", "CTACC", "CTACG", "CTACT", "CTAGA", "CTAGC", "CTAGG", "CTAGT", "CTATA", "CTATC", "CTATG", "CTATT", "CTCAA", "CTCAC", "CTCAG", "CTCAT", "CTCCA", "CTCCC", "CTCCG", "CTCCT", "CTCGA", "CTCGC", "CTCGG", "CTCGT", "CTCTA", "CTCTC", "CTCTG", "CTCTT", "CTGAA", "CTGAC", "CTGAG", "CTGAT", "CTGCA", "CTGCC", "CTGCG", "CTGCT", "CTGGA", "CTGGC", "CTGGG", "CTGGT", "CTGTA", "CTGTC", "CTGTG", "CTGTT", "CTTAA", "CTTAC", "CTTAG", "CTTAT", "CTTCA", "CTTCC", "CTTCG", "CTTCT", "CTTGA", "CTTGC", "CTTGG", "CTTGT", "CTTTA", "CTTTC", "CTTTG", "CTTTT"}},
								{'G', { "G", "GA", "GC", "GG", "GT", "GAA", "GAC", "GAG", "GAT", "GCA", "GCC", "GCG", "GCT", "GGA", "GGC", "GGG", "GGT", "GTA", "GTC", "GTG", "GTT", "GAAA", "GAAC", "GAAG", "GAAT", "GACA", "GACC", "GACG", "GACT", "GAGA", "GAGC", "GAGG", "GAGT", "GATA", "GATC", "GATG", "GATT", "GCAA", "GCAC", "GCAG", "GCAT", "GCCA", "GCCC", "GCCG", "GCCT", "GCGA", "GCGC", "GCGG", "GCGT", "GCTA", "GCTC", "GCTG", "GCTT", "GGAA", "GGAC", "GGAG", "GGAT", "GGCA", "GGCC", "GGCG", "GGCT", "GGGA", "GGGC", "GGGG", "GGGT", "GGTA", "GGTC", "GGTG", "GGTT", "GTAA", "GTAC", "GTAG", "GTAT", "GTCA", "GTCC", "GTCG", "GTCT", "GTGA", "GTGC", "GTGG", "GTGT", "GTTA", "GTTC", "GTTG", "GTTT", "GAAAA", "GAAAC", "GAAAG", "GAAAT", "GAACA", "GAACC", "GAACG", "GAACT", "GAAGA", "GAAGC", "GAAGG", "GAAGT", "GAATA", "GAATC", "GAATG", "GAATT", "GACAA", "GACAC", "GACAG", "GACAT", "GACCA", "GACCC", "GACCG", "GACCT", "GACGA", "GACGC", "GACGG", "GACGT", "GACTA", "GACTC", "GACTG", "GACTT", "GAGAA", "GAGAC", "GAGAG", "GAGAT", "GAGCA", "GAGCC", "GAGCG", "GAGCT", "GAGGA", "GAGGC", "GAGGG", "GAGGT", "GAGTA", "GAGTC", "GAGTG", "GAGTT", "GATAA", "GATAC", "GATAG", "GATAT", "GATCA", "GATCC", "GATCG", "GATCT", "GATGA", "GATGC", "GATGG", "GATGT", "GATTA", "GATTC", "GATTG", "GATTT", "GCAAA", "GCAAC", "GCAAG", "GCAAT", "GCACA", "GCACC", "GCACG", "GCACT", "GCAGA", "GCAGC", "GCAGG", "GCAGT", "GCATA", "GCATC", "GCATG", "GCATT", "GCCAA", "GCCAC", "GCCAG", "GCCAT", "GCCCA", "GCCCC", "GCCCG", "GCCCT", "GCCGA", "GCCGC", "GCCGG", "GCCGT", "GCCTA", "GCCTC", "GCCTG", "GCCTT", "GCGAA", "GCGAC", "GCGAG", "GCGAT", "GCGCA", "GCGCC", "GCGCG", "GCGCT", "GCGGA", "GCGGC", "GCGGG", "GCGGT", "GCGTA", "GCGTC", "GCGTG", "GCGTT", "GCTAA", "GCTAC", "GCTAG", "GCTAT", "GCTCA", "GCTCC", "GCTCG", "GCTCT", "GCTGA", "GCTGC", "GCTGG", "GCTGT", "GCTTA", "GCTTC", "GCTTG", "GCTTT", "GGAAA", "GGAAC", "GGAAG", "GGAAT", "GGACA", "GGACC", "GGACG", "GGACT", "GGAGA", "GGAGC", "GGAGG", "GGAGT", "GGATA", "GGATC", "GGATG", "GGATT", "GGCAA", "GGCAC", "GGCAG", "GGCAT", "GGCCA", "GGCCC", "GGCCG", "GGCCT", "GGCGA", "GGCGC", "GGCGG", "GGCGT", "GGCTA", "GGCTC", "GGCTG", "GGCTT", "GGGAA", "GGGAC", "GGGAG", "GGGAT", "GGGCA", "GGGCC", "GGGCG", "GGGCT", "GGGGA", "GGGGC", "GGGGG", "GGGGT", "GGGTA", "GGGTC", "GGGTG", "GGGTT", "GGTAA", "GGTAC", "GGTAG", "GGTAT", "GGTCA", "GGTCC", "GGTCG", "GGTCT", "GGTGA", "GGTGC", "GGTGG", "GGTGT", "GGTTA", "GGTTC", "GGTTG", "GGTTT", "GTAAA", "GTAAC", "GTAAG", "GTAAT", "GTACA", "GTACC", "GTACG", "GTACT", "GTAGA", "GTAGC", "GTAGG", "GTAGT", "GTATA", "GTATC", "GTATG", "GTATT", "GTCAA", "GTCAC", "GTCAG", "GTCAT", "GTCCA", "GTCCC", "GTCCG", "GTCCT", "GTCGA", "GTCGC", "GTCGG", "GTCGT", "GTCTA", "GTCTC", "GTCTG", "GTCTT", "GTGAA", "GTGAC", "GTGAG", "GTGAT", "GTGCA", "GTGCC", "GTGCG", "GTGCT", "GTGGA", "GTGGC", "GTGGG", "GTGGT", "GTGTA", "GTGTC", "GTGTG", "GTGTT", "GTTAA", "GTTAC", "GTTAG", "GTTAT", "GTTCA", "GTTCC", "GTTCG", "GTTCT", "GTTGA", "GTTGC", "GTTGG", "GTTGT", "GTTTA", "GTTTC", "GTTTG", "GTTTT"}},
								{'T', { "T", "TA", "TC", "TG", "TT", "TAA", "TAC", "TAG", "TAT", "TCA", "TCC", "TCG", "TCT", "TGA", "TGC", "TGG", "TGT", "TTA", "TTC", "TTG", "TTT", "TAAA", "TAAC", "TAAG", "TAAT", "TACA", "TACC", "TACG", "TACT", "TAGA", "TAGC", "TAGG", "TAGT", "TATA", "TATC", "TATG", "TATT", "TCAA", "TCAC", "TCAG", "TCAT", "TCCA", "TCCC", "TCCG", "TCCT", "TCGA", "TCGC", "TCGG", "TCGT", "TCTA", "TCTC", "TCTG", "TCTT", "TGAA", "TGAC", "TGAG", "TGAT", "TGCA", "TGCC", "TGCG", "TGCT", "TGGA", "TGGC", "TGGG", "TGGT", "TGTA", "TGTC", "TGTG", "TGTT", "TTAA", "TTAC", "TTAG", "TTAT", "TTCA", "TTCC", "TTCG", "TTCT", "TTGA", "TTGC", "TTGG", "TTGT", "TTTA", "TTTC", "TTTG", "TTTT", "TAAAA", "TAAAC", "TAAAG", "TAAAT", "TAACA", "TAACC", "TAACG", "TAACT", "TAAGA", "TAAGC", "TAAGG", "TAAGT", "TAATA", "TAATC", "TAATG", "TAATT", "TACAA", "TACAC", "TACAG", "TACAT", "TACCA", "TACCC", "TACCG", "TACCT", "TACGA", "TACGC", "TACGG", "TACGT", "TACTA", "TACTC", "TACTG", "TACTT", "TAGAA", "TAGAC", "TAGAG", "TAGAT", "TAGCA", "TAGCC", "TAGCG", "TAGCT", "TAGGA", "TAGGC", "TAGGG", "TAGGT", "TAGTA", "TAGTC", "TAGTG", "TAGTT", "TATAA", "TATAC", "TATAG", "TATAT", "TATCA", "TATCC", "TATCG", "TATCT", "TATGA", "TATGC", "TATGG", "TATGT", "TATTA", "TATTC", "TATTG", "TATTT", "TCAAA", "TCAAC", "TCAAG", "TCAAT", "TCACA", "TCACC", "TCACG", "TCACT", "TCAGA", "TCAGC", "TCAGG", "TCAGT", "TCATA", "TCATC", "TCATG", "TCATT", "TCCAA", "TCCAC", "TCCAG", "TCCAT", "TCCCA", "TCCCC", "TCCCG", "TCCCT", "TCCGA", "TCCGC", "TCCGG", "TCCGT", "TCCTA", "TCCTC", "TCCTG", "TCCTT", "TCGAA", "TCGAC", "TCGAG", "TCGAT", "TCGCA", "TCGCC", "TCGCG", "TCGCT", "TCGGA", "TCGGC", "TCGGG", "TCGGT", "TCGTA", "TCGTC", "TCGTG", "TCGTT", "TCTAA", "TCTAC", "TCTAG", "TCTAT", "TCTCA", "TCTCC", "TCTCG", "TCTCT", "TCTGA", "TCTGC", "TCTGG", "TCTGT", "TCTTA", "TCTTC", "TCTTG", "TCTTT", "TGAAA", "TGAAC", "TGAAG", "TGAAT", "TGACA", "TGACC", "TGACG", "TGACT", "TGAGA", "TGAGC", "TGAGG", "TGAGT", "TGATA", "TGATC", "TGATG", "TGATT", "TGCAA", "TGCAC", "TGCAG", "TGCAT", "TGCCA", "TGCCC", "TGCCG", "TGCCT", "TGCGA", "TGCGC", "TGCGG", "TGCGT", "TGCTA", "TGCTC", "TGCTG", "TGCTT", "TGGAA", "TGGAC", "TGGAG", "TGGAT", "TGGCA", "TGGCC", "TGGCG", "TGGCT", "TGGGA", "TGGGC", "TGGGG", "TGGGT", "TGGTA", "TGGTC", "TGGTG", "TGGTT", "TGTAA", "TGTAC", "TGTAG", "TGTAT", "TGTCA", "TGTCC", "TGTCG", "TGTCT", "TGTGA", "TGTGC", "TGTGG", "TGTGT", "TGTTA", "TGTTC", "TGTTG", "TGTTT", "TTAAA", "TTAAC", "TTAAG", "TTAAT", "TTACA", "TTACC", "TTACG", "TTACT", "TTAGA", "TTAGC", "TTAGG", "TTAGT", "TTATA", "TTATC", "TTATG", "TTATT", "TTCAA", "TTCAC", "TTCAG", "TTCAT", "TTCCA", "TTCCC", "TTCCG", "TTCCT", "TTCGA", "TTCGC", "TTCGG", "TTCGT", "TTCTA", "TTCTC", "TTCTG", "TTCTT", "TTGAA", "TTGAC", "TTGAG", "TTGAT", "TTGCA", "TTGCC", "TTGCG", "TTGCT", "TTGGA", "TTGGC", "TTGGG", "TTGGT", "TTGTA", "TTGTC", "TTGTG", "TTGTT", "TTTAA", "TTTAC", "TTTAG", "TTTAT", "TTTCA", "TTTCC", "TTTCG", "TTTCT", "TTTGA", "TTTGC", "TTTGG", "TTTGT", "TTTTA", "TTTTC", "TTTTG", "TTTTT"}}
};

/* Checks that the filepath is readable and exits if it is not. */
static inline void assert_readable(const string& path) {
	if (access(path.c_str(), R_OK) == -1) {
		std::cerr << PROGRAM ": error: `" << path << "': " << strerror(errno) << std::endl; 
		exit(EXIT_FAILURE); 
	}
}

/* Checks if the base is ATGC. */
bool isAcceptedBase(unsigned char C) {
	return (C=='A' || C=='T' || C=='G' || C=='C'); 
}

/* Find the first only ATGC kmer starting at the beginning of a sequence.
 * 	Assumption: no insertions or deletions. */
unsigned findFirstAcceptedKmer(unsigned b_i, const string& contigSeq) {
	for (unsigned i=b_i; i+opt::k<contigSeq.size();) {
		if (isAcceptedBase(toupper(contigSeq.at(i)))) {
			bool good_kmer=true; 
			for (unsigned j=i+1; j<i+opt::k; j++) {
				if (!isAcceptedBase(toupper(contigSeq.at(j)))) {
					good_kmer=false;
					i=j+1; 
					break;
				}
			}
			if (good_kmer) {
				return i; 
			}
		} else {
			i++; 
		}
	}
	return contigSeq.size()-1; 
}

/* Helper for filling out the LPS array for detecting a low complexity repeat. */
void computeLPSArray(string possible_repeat, int n, vector<int>& lps) {
	int len=0; 
	int i; 

	lps[0] = 0; 
	i=1;
	while (i<n) {
		if (possible_repeat[i] == possible_repeat[len]) {
			len++; 
			lps[i] = len; 
			i++;
		} else {
			if (len!=0)
				len = lps[len-1]; 
			else {
				lps[i] = 0; 
				i++;
			}
		}
	}
}

/* Determines if a string is a low complexity repeat of a word. */
bool isRepeatInsertion(string possible_repeat) {
	int n = possible_repeat.size(); 
	vector<int> lps(n); 

	computeLPSArray(possible_repeat, n, lps); 

	int len = lps[n-1]; 
	return (len>0 && n%(n-len) == 0);
}

/* Struct that keeps track of details for substitutions. */
struct sRec {
	unsigned pos; 
	unsigned char draft_char;
	unsigned char sub_base; 
	unsigned num_support;
};

/* START: RopeLink Structure and Functions -- todo: abstract this into its own class. */
struct RopeLink {
	RopeLink *left, *right; // left== refers to parent; right== refers to following
	int node_type=-1; // -1==unset; 0==position; 1==character;
	size_t s_pos; // start position of a position node
	size_t e_pos; // end position of a position node
	unsigned char c; // character for a character node
	unsigned num_support=0; // check_present for an insertion/deletion
};

/* Creates a position RopeLink at <node> with parent <left> and s_pos <s_pos> and e_pos <e_pos>. */
void createPositionRopeLink(RopeLink *& node, RopeLink *& left, size_t s_pos, size_t e_pos) {
	RopeLink * tmp = (struct RopeLink* ) malloc(sizeof(struct RopeLink)); 
	tmp->node_type = 0; 
	tmp->left = left;
	if (left != nullptr) left->right = tmp;
	tmp->right = nullptr;
	tmp->s_pos = s_pos; 
	tmp->e_pos = e_pos;

	node = tmp; 
}

/* Creates a character RopeLink at <node> with parent <left> and character <c>. */
void createCharacterRopeLink(RopeLink *& node, RopeLink *& left, unsigned char c) {
	RopeLink *tmp = (struct RopeLink* ) malloc(sizeof(struct RopeLink)); 
	tmp->node_type = 1;
	tmp->left = left; 
	if (left != nullptr) left->right = tmp;
	tmp->right = nullptr; 
	tmp->c = c; 

	node = tmp; 
}

/* Removes a RopeLink and adjusts the parent and next node pointers. */
void removeRopeLink(RopeLink *& to_remove) {
	RopeLink * prev_node = to_remove->left; 
	if (prev_node != nullptr) {
		prev_node->right = to_remove->right; 
	}
	RopeLink * next_node = to_remove->right; 
	if (next_node != nullptr) {
		next_node->left = prev_node; 
	}
	to_remove->left = nullptr; 
	to_remove->right = nullptr;
	to_remove->node_type = -1; 	

	free(to_remove); 
}

/* Makes a character insertion RIGHT BEFORE <insert_pos> by creating a character node holding <c> with <num_support> 
 * 	- sets <node> to your insertion node (the node that holds the character <c>). */
void makeCharacterInsertion(RopeLink *& node, int insert_pos, unsigned char c, unsigned num_support) {
	RopeLink *insertion_node = nullptr; 
	RopeLink *prev_node = nullptr; 
	RopeLink *next_node = nullptr; 
	if (node->node_type == 0) {
		if (insert_pos <= node->s_pos) {
			// insertion occurs at the beginning of this node 
			// 	which means we've made an insertion/deletion or edit at this area
			//		BEFORE: .... <------> original position node <-------> .....
			//		AFTER: ..... <----> insertion_node <--------> original position node <------> ....
			prev_node = node->left; 
			createCharacterRopeLink(insertion_node, prev_node, c); 
			next_node = node; 
			insertion_node->right = next_node; 
			if (next_node != nullptr) next_node->left = insertion_node; 
		} else {
			// insertion occurs after the beginning of this node and in the middle of the position node
			// 	BEFORE: ......<---------> original position node <------> .....
			// 	AFTER:  ..... <--> part 1 of position node <--> insertion_node <--> part 2 of position node <---> ...
			size_t tmp_epos = node->e_pos; 
			RopeLink *tmp_right_flank = node->right;
			prev_node = node;
			prev_node->e_pos = insert_pos-1; 
			createCharacterRopeLink(insertion_node, prev_node, c); 
			createPositionRopeLink(next_node, insertion_node, insert_pos, tmp_epos); 
			next_node->right = tmp_right_flank; 
			if (tmp_right_flank != nullptr) tmp_right_flank->left = next_node;
		}
	} else if (node->node_type == 1) {
		// make a character insertion before a character
		// 	BEFORE: .... <-----> original character node <----> ....
		// 	AFTER: .... <------> insertion_node <----> original character node <----> .....
		prev_node = node->left; 
		next_node = node;
		createCharacterRopeLink(insertion_node, prev_node, c); 
		insertion_node->right = next_node; 
		if (next_node != nullptr) next_node->left = insertion_node; 
	}
	insertion_node->num_support = num_support; 
	// adjust our tNode
	node = insertion_node;
}

/* Makes a deletion by deleting the characters starting and including <delete_pos> of size <num_del>. */
void makeDeletion(RopeLink *&node, unsigned& delete_pos, unsigned num_del, unsigned num_support) {
	RopeLink * tmp_right_node = node->right;
	// making a deletion from a position node
	if (node->node_type == 0) {
		// save the old end position and set it to the new one
		size_t tmp_epos = node->e_pos; 
		unsigned leftover_del=0; 
		// deletion occurs at the beginning of this node
		if (delete_pos <= node->s_pos) {
			if (delete_pos+num_del <= node->e_pos)  {
				// deletion occurs at the beginning of the node and is kept within the node
				node->s_pos = delete_pos+num_del; 
				node->num_support = num_support;
				delete_pos += num_del;
				return;
			} else {
				// deletion occurs at the beginning of the node and goes past the current node
				leftover_del = delete_pos+num_del-node->e_pos;
				// deletion removed this entire node
				RopeLink *to_remove = node; 
				node = node->right; 
				removeRopeLink(to_remove); 
			}
		} else {
			if (delete_pos+num_del <= node->e_pos) {
				// deletion occurs starting in the middle of the node and is kept within the node
				node->e_pos = delete_pos-1; 
				RopeLink *split_node = nullptr; 
				createPositionRopeLink(split_node, node, delete_pos+num_del, tmp_epos); 
				split_node->num_support = num_support;
				node = node->right; 
				delete_pos += num_del; 
				return;
			} else {
				// deletion occurs starting in the middle of the node and continues past the node
				leftover_del = delete_pos+num_del-node->e_pos; 
				delete_pos = node->e_pos+1; 
				node->e_pos = delete_pos-1; 
				node = node->right; 
			}
		}
		// if there are still deletions to finish continue because it has continued past this node
		if (leftover_del > 0) {
			if (node != nullptr) {
				if (node->node_type == 0) delete_pos = node->s_pos; 
				else if (node->node_type == 1) delete_pos = tmp_epos++; 
				makeDeletion(node, delete_pos, leftover_del, num_support); 
			}
		}
	// we're making a deletion from a character node
	} else if (node->node_type == 1) {
		// delete this node
		RopeLink *to_remove = node; 
		node = node->right; 
		removeRopeLink(to_remove); 
		num_del--; 
		if (node != nullptr && num_del > 0) {
			if (node->node_type == 0) {
				delete_pos = node->s_pos; 
				makeDeletion(node, delete_pos, num_del, num_support); 
			} else if (node->node_type == 1) {
				makeDeletion(node, delete_pos, num_del, num_support); 
			}
		} 
	}
}

/* Starts at root and removes all of the nodes following it. */
void cleanRopeLinks(RopeLink *&master_root) {
	RopeLink * to_remove;
	do {
		to_remove = master_root; 
		master_root = master_root->right; 
		removeRopeLink(to_remove); 
	} while (master_root != nullptr); 
}
/* END: RopeLink Struct and Functions */

/* Returns the character at pos based on the RopeLink structure. */
unsigned char getCharacter(unsigned pos, RopeLink *node, const string& contigSeq) {
	if (node->node_type == 0) return contigSeq.at(pos); 
	else if (node->node_type == 1) return node->c; 
	unsigned char c;
	return c; 
}

/* Increments the position and adjusts the node accordingly based on RopeLink structure. */
void increment(unsigned& pos, RopeLink *& node) {
	if (node->node_type == 0) {
		pos++; 
		if (pos > node->e_pos) {
			node = node->right; 
			if (node != nullptr && node->node_type == 0 && pos < node->s_pos)
				pos = node->s_pos;
		} 
	} else if (node->node_type == 1) {
		node = node->right; 
		if (node != nullptr && node->node_type == 0 && pos < node->s_pos)
			pos = node->s_pos;
	}
}

/* Find the first accepted kmer (contains only ATGC characters) starting anywhere, based on the RopeLink structure. */
string findAcceptedKmer(unsigned& h_seq_i, unsigned& t_seq_i,
			RopeLink*& hNode, RopeLink *& tNode, const string& contigSeq) {
	// temporary values
	string kmer_str; 	
	RopeLink * curr_node = tNode;
	RopeLink * new_hNode, * new_tNode;
	unsigned i=t_seq_i; 
	while (i<contigSeq.size() && curr_node != nullptr) {
		unsigned char c; 
		c = getCharacter(i, curr_node, contigSeq); 
		if (isAcceptedBase(toupper(c))) {
			string kmer_str; 
			kmer_str += c; 
			std::cout << "starting kmer: " << kmer_str << " " << i << std::endl; 
			new_hNode = curr_node; 
			unsigned j=i; 
			increment(j, curr_node); 
			// continue until you cant
			while (j<contigSeq.size() && curr_node != nullptr) {
				c = getCharacter(j, curr_node, contigSeq); 
				if (!isAcceptedBase(toupper(c))) {
					i=j; 
					break;
				} 
				kmer_str += c;
				std::cout << kmer_str << std::endl; 
			        if (kmer_str.size() == opt::k) break;	
				increment(j, curr_node); 
			}
			std::cout << "finished looking for kmer starting here: " << kmer_str << std::endl; 
			// you found a good kmer so return it and adjust
			if (kmer_str.size() == opt::k) {
				h_seq_i = i; 
				t_seq_i = j; 
				hNode = new_hNode;
				tNode = curr_node;
				std::cout << "set everything: " << h_seq_i << " " << t_seq_i << " " << kmer_str << std::endl; 
				return kmer_str; 				
			}
		}
		increment(i, curr_node); 
	}

	h_seq_i = contigSeq.length(); 
	t_seq_i = contigSeq.length(); 
	return ""; 
}

/* Get the previous insertion (aka continuous string of character nodes) starting at this node tNode. */
string getPrevInsertion(unsigned t_seq_i, RopeLink *tNode) {
	// if we just finished the insertion
	if (tNode != nullptr && tNode->node_type == 0 && t_seq_i == tNode->s_pos) {
		tNode = tNode->left; 
	}
	string prev_insertion;
	while (tNode != nullptr && tNode->node_type == 1) {
		prev_insertion += tNode->c; 
		tNode = tNode->left;
	}
	return prev_insertion;
}

/* Write the edits and new draft contig into respective files. */
void writeEditsToFile(FILE* dfout, FILE* rfout, const string& contigHdr, const string& contigSeq,
		RopeLink *& root, queue<sRec>& substitution_record) {
	fprintf(dfout, ">%s\n", contigHdr.c_str()); 
	RopeLink * curr_node = root;
	string insertion_bases=""; 
	int num_support=-1; 
	unsigned char draft_char;
	unsigned pos; 
	// track a deletion
	string deleted_bases="";
	while (curr_node != nullptr) {
		if (curr_node->node_type == 0) {
			draft_char = contigSeq.at(curr_node->s_pos);
			// log an insertion if it occured before this
			if (!insertion_bases.empty()) {
				fprintf(rfout, "%s\t%d\t%c\t%c%s\t%d\n",
						contigHdr.c_str(), pos+1, draft_char,'+', insertion_bases.c_str(), num_support); 
				insertion_bases = "";  
				num_support = -1; 
			}
			// log all the substitutions up to this point
			while (!substitution_record.empty() && substitution_record.front().pos <= curr_node->e_pos){
				fprintf(rfout, "%s\t%d\t%c\t%c\t%d\n",
						contigHdr.c_str(), substitution_record.front().pos+1,
						substitution_record.front().draft_char,
						substitution_record.front().sub_base,
						substitution_record.front().num_support);
				substitution_record.pop();
			} 
			fprintf(dfout, "%s", contigSeq.substr(curr_node->s_pos, (curr_node->e_pos-curr_node->s_pos+1)).c_str()); 
			pos = curr_node->e_pos+1; 
		} else if (curr_node->node_type == 1) {
			insertion_bases += curr_node->c; 
			if (num_support==-1) num_support = curr_node->num_support; 
			fprintf(dfout, "%c", curr_node->c); 
		}
		curr_node = curr_node->right; 
		if (curr_node != nullptr && curr_node->node_type == 0 && curr_node->s_pos != pos) {
			// print out the deletion
			fprintf(rfout, "%s\t%d\t%c\t%c%s\t%d\n", 
					contigHdr.c_str(), pos+1, contigSeq.at(pos), '-', 
					contigSeq.substr(pos, (curr_node->s_pos-pos)).c_str(), curr_node->num_support); 
		}
	}
	fprintf(dfout, "\n");
}

/* Roll ntHash using the Rope Data structure. */
bool roll(unsigned& h_seq_i, unsigned& t_seq_i, 
		RopeLink *& hNode, RopeLink *&tNode,
		uint64_t& fhVal, uint64_t& rhVal, uint64_t *&hVal,
		const string& contigSeq, unsigned char& charIn) {

	// quit if h_seq_i is out of scope
	if (h_seq_i >= contigSeq.size() || hNode == nullptr) return false;

	unsigned char charOut = getCharacter(h_seq_i, hNode, contigSeq); 
	increment(h_seq_i, hNode); 
	if (hNode == nullptr) return false; 

	increment(t_seq_i, tNode); 
	if (t_seq_i >= contigSeq.size() || tNode == nullptr) return false; 
	charIn = getCharacter(t_seq_i, tNode, contigSeq); 

	// roll the hash
	NTMC64(charOut, charIn, opt::k, opt::h, fhVal, rhVal, hVal); 
	return true;
}

/* Returns a substring using the RopeLink structure starting at <start_pos> of length <length>. */
string getSubstring(unsigned start_pos, RopeLink *sNode, int length, const string& contigSeq) {
	string substring; 
	substring.reserve(length);
	stringstream ss(substring); 
	unsigned i=start_pos;
	int len = length; 
	while (len > 0 && i < contigSeq.size() &&  sNode != nullptr) {
		if (sNode->node_type == 0) {
			if (i < sNode->s_pos) i = sNode->s_pos;
			ss << contigSeq.at(i); 
			len--;
			i++; 
			if (i > sNode->e_pos) sNode = sNode->right; 
		} else if (sNode->node_type == 1) {
			ss << sNode->c;
			len--;	
			sNode = sNode->right;
		}
	}
	return ss.str(); 
}

/* Try a deletion in ntEdit. */
int tryDeletion(const unsigned char draft_char, unsigned num_deletions,
		unsigned& h_seq_i, unsigned& t_seq_i, 
		RopeLink *&hNode, RopeLink *&tNode,
		uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal, 
		string& before_t, string& after_t, BloomFilter& bloom) {

	// set temporary values
	uint64_t temp_fhVal = fhVal, temp_rhVal = rhVal;
	unsigned temp_h_seq_i = 0, temp_t_seq_i = opt::k-1; 
	RopeLink *temp_hNode = hNode, *temp_tNode = tNode; 
	string deleted_bases; 

	// make the subsetstring we are iterating on
	string deletion_substring; 
	deletion_substring.reserve((2*opt::k)-1-num_deletions); 
	stringstream ss(deletion_substring); 
	ss << before_t << after_t.substr(num_deletions); 
	deletion_substring = ss.str();
	deleted_bases = after_t.substr(0, num_deletions); 
	NTMC64_changelast(draft_char, deletion_substring.at(temp_t_seq_i), opt::k, opt::h, temp_fhVal, temp_rhVal, hVal); 

	// verify the deletion with a subset
	unsigned check_present=0; 
	if (bloom.contains(hVal)) check_present++; // check for changing the kmer after deletion
	for (unsigned k=1; k<=(opt::k-2) && temp_h_seq_i<deletion_substring.size() 
			&& ++temp_t_seq_i <deletion_substring.size(); k++) {
		NTMC64(deletion_substring.at(temp_h_seq_i), deletion_substring.at(temp_t_seq_i), opt::k, opt::h,
				temp_fhVal, temp_rhVal, hVal); 
		temp_h_seq_i++; 
		if (k%3 == 0 && bloom.contains(hVal)) check_present++;
	}
	
	if (opt::verbose)
		std::cout << "\t\tdeleting: " << deleted_bases << " check_present: " << check_present << std::endl; 
	if (check_present >= ((float) opt::k / opt::edit_threshold)) {
		return check_present;
	} 
	return 0;
}

/* Try indel combinations starting with index_char. */
bool tryIndels(const unsigned char draft_char, const unsigned char index_char, 
		unsigned char& charIn,
		unsigned& num_deletions,
		unsigned& h_seq_i, unsigned& t_seq_i,
		RopeLink *&hNode, RopeLink *& tNode,
		uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal,
		string& before_t, string& after_t, const string& contigSeq, 
		BloomFilter& bloom, unsigned& best_edit_type, string& best_indel, unsigned& best_num_support) {

	// initialize temporary values
	uint64_t temp_fhVal, temp_rhVal; 
	unsigned temp_h_seq_i, temp_t_seq_i; 
	unsigned temp_best_num_support=0;
	string temp_best_indel; 
	unsigned temp_best_edit_type=0;


	// try all of the combinations of indels starting with our index_char
	for (unsigned i=0; i<num_tries[opt::max_insertions]; i++) {
		// gather the insertion bases
		string insertion_bases = multi_possible_bases[index_char][i]; 

		// subset string
		string subset_str; 
		subset_str.reserve((2*opt::k)-1+insertion_bases.size()); 
		stringstream ss(subset_str); 
		ss << before_t << insertion_bases << after_t; 
		subset_str = ss.str();
		//std::cout << subset_str << std::endl; 

		// set temporary values
		temp_fhVal = fhVal;
		temp_rhVal = rhVal;
		temp_h_seq_i = 0; 
		temp_t_seq_i = opt::k-1; 
		
		// change the last base
		NTMC64_changelast(draft_char, index_char, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal); 
		unsigned check_present=0; 
		for (unsigned k=1; k<=opt::k-1 && temp_h_seq_i<subset_str.size()
				&& ++temp_t_seq_i<subset_str.size(); k++) {
			NTMC64(subset_str.at(temp_h_seq_i), subset_str.at(temp_t_seq_i), opt::k, opt::h, 
					temp_fhVal, temp_rhVal, hVal);
			temp_h_seq_i++; 
			if (k%3 == 1 && bloom.contains(hVal))
				check_present++; 
		}
		if (opt::verbose)
			std::cout << "\t\tinserting: " << insertion_bases << " check_present: " << check_present << std::endl;
		// if the insertion is good, store the insertion accordingly 
		if (check_present >= ((float) opt::k / opt::edit_threshold)) {
			if (opt::mode == 0) {
				// if we are in default mode, we just accept this first good insertion and return
				best_edit_type = 2; 
				best_indel = insertion_bases; 
				best_num_support = check_present; 
				return true;
			} else if (opt::mode == 1 || opt::mode == 2) {
				// if we are in some deep mode, we look for the best indel within index char first
				if (check_present > temp_best_num_support) {
					temp_best_edit_type = 2; 
					temp_best_indel = insertion_bases; 
					temp_best_num_support = check_present; 
				}

			}
		}

		if (num_deletions <= opt::max_deletions) {
			unsigned del_support = tryDeletion(draft_char, num_deletions, h_seq_i, t_seq_i, 
								hNode, tNode, fhVal, rhVal, hVal, before_t, after_t, bloom);
			if (del_support > 0) {
				if (opt::mode == 0) {
					best_edit_type = 3;
					best_indel = after_t.substr(0, num_deletions); 
					best_num_support = del_support;
					return true;
				} else if (opt::mode == 1 || opt::mode == 2) {
					if (del_support > temp_best_num_support) {
						temp_best_edit_type = 3;
						temp_best_indel = after_t.substr(0, num_deletions);
						temp_best_num_support = del_support;
					}
				}
			}
			num_deletions++; 
		}
	}
	
	// report the best indel info
	if (temp_best_num_support > 0) {
		if ((opt::mode == 2 && temp_best_num_support > best_num_support) || opt::mode == 1) {
			best_edit_type = temp_best_edit_type;
			best_indel = temp_best_indel; 
			best_num_support = temp_best_num_support;
		} 
		return true;
	}
	return false; 
}

/* Kmerize and polish the contig. */
void kmerizeAndCorrect(string& contigHdr, string& contigSeq, unsigned seqLen, BloomFilter& bloom, 
		FILE* dfout, FILE* rfout) {

	// initialize values for hashing
	uint64_t fhVal, rhVal;
	uint64_t* hVal;
	unsigned char charIn, draft_char; 
	hVal = new uint64_t[opt::h]; 


	// vector to record substitutions
	queue<sRec> substitution_record; 
	
	// initialize and readjust the first character depending on the first N or nonATGC kmer
	unsigned h_seq_i = findFirstAcceptedKmer(0, contigSeq); 
	unsigned t_seq_i = h_seq_i+opt::k-1;

	// intialize our seed kmer
	if (h_seq_i+opt::k-1 < seqLen) {
		NTMC64(contigSeq.substr(h_seq_i, opt::k).c_str(), opt::k, opt::h, fhVal, rhVal, hVal); 
		charIn=contigSeq.at(t_seq_i); 
	} 	
	
	// initialize our master root
	RopeLink * master_root = nullptr;
	RopeLink * null_parent = nullptr;
	createPositionRopeLink(master_root, null_parent, 0, seqLen-1); 
	// pointers to nodes for h_seq_i and t_seq_i
	RopeLink * hNode = master_root;
	RopeLink * tNode = master_root;

	bool continue_edit = true;
	do {
		if (h_seq_i+opt::k-1 >= seqLen) break;
		if (opt::verbose) 
			std::cout << h_seq_i << " " << t_seq_i << " " << charIn << " " << hVal[0] << hVal[1] << hVal[2] << std::endl;

		if (!bloom.contains(hVal)) {
			// make temporary value holders
			uint64_t temp_fhVal = fhVal; 
			uint64_t temp_rhVal = rhVal; 
			unsigned temp_h_seq_i = h_seq_i; 
			unsigned temp_t_seq_i = t_seq_i; 
			RopeLink* temp_hNode = hNode;
			RopeLink* temp_tNode = tNode;

			// set draft char if we choose to make an edit later
			draft_char = toupper(charIn); 

			// confirm missing by checking subset
			unsigned check_missing=0; 
			bool do_not_fix = false; 
			for (unsigned k=1; k<=opt::k && temp_h_seq_i<seqLen; k++) {
				if (roll(temp_h_seq_i, temp_t_seq_i, temp_hNode, temp_tNode,
							temp_fhVal, temp_rhVal, hVal, contigSeq, charIn)) {
					if (!isAcceptedBase(toupper(charIn))) {
						do_not_fix = true; 
						break;
					}
					if (k%3 == 1 && !bloom.contains(hVal)) check_missing++; 
				} else break;
			}

			if (opt::verbose) 
				std::cout << "\tcheck_missing: " << check_missing << std::endl; 
			if (!do_not_fix && check_missing>=((float) opt::k / opt::missing_threshold)) {
				// recorders
				unsigned num_deletions=1;
				unsigned best_edit_type=0; // 0=no edit made; 1=substitution; 2=insertion; 3=deletions
				string best_indel; 
				unsigned char best_sub_base; 
				unsigned best_num_support=0; 
				
				// everything before t_seq_i
				string before_t = getSubstring(h_seq_i, hNode, opt::k-1, contigSeq); 
				// everything including t_seq_i and after
				string after_t = getSubstring(t_seq_i, tNode, opt::k+opt::max_deletions, contigSeq); 

				// try substitution 
				for (const unsigned char sub_base : bases_array[draft_char]) {
					//reset the temporary values
					temp_fhVal = fhVal; 
					temp_rhVal = rhVal; 
					temp_hNode = hNode; 
					temp_tNode = tNode;

					// hash the substitution change
					NTMC64_changelast(draft_char, sub_base, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal); 

					// only do verification of substition if it is found in bloom filter
					if (bloom.contains(hVal) || opt::mode == 2) {
						// reset temporary values 
						temp_h_seq_i = h_seq_i; 
						temp_t_seq_i = t_seq_i; 

						// change the substitution 
						if (temp_tNode->node_type == 0) contigSeq.at(temp_t_seq_i) = sub_base; 
						else if (temp_tNode->node_type == 1) temp_tNode->c = sub_base;

						// check the subset to see if this substition is good
						unsigned check_present=0; 
						for (unsigned k=1; k<=opt::k && temp_h_seq_i < seqLen 
								&& temp_t_seq_i<seqLen; k++) {
							if (roll(temp_h_seq_i, temp_t_seq_i, temp_hNode, temp_tNode,
									temp_fhVal, temp_rhVal, hVal, contigSeq, charIn)) {
								if (k%3 == 1 && bloom.contains(hVal)) check_present++; 
							} else break;
						}

						// revert the substitution
						if (tNode->node_type == 0) contigSeq.at(t_seq_i) = draft_char; 
						else if (tNode->node_type == 1) tNode->c = draft_char; 

						if (opt::verbose)
							std::cout << "\t\tsub: " << sub_base << " check_present: " 
								<< check_present << std::endl; 
						if (check_present >= ((float) opt::k / opt::edit_threshold)) {
							// update the best substitution
							if (check_present > best_num_support) {
								best_edit_type = 1;
								best_sub_base = sub_base; 
								best_num_support = check_present; 
							}
							// if we aren't exhaustively trying all edit combinations, 
							// 	then just do substitutions from now on
							if (opt::mode == 0 || opt::mode == 1) {
								continue;
							}
						}
						// if we are exhaustively trying all edit combinations or 
						// 	havent found a good substitution yet, then try indels
						if (opt::mode == 2 || best_edit_type != 1) {
							if (tryIndels(draft_char, sub_base, charIn, num_deletions, h_seq_i, t_seq_i,
									      	hNode, tNode,
									fhVal, rhVal, hVal, before_t, after_t, contigSeq, bloom,
									best_edit_type, best_indel, best_num_support)) {
								if (opt::mode == 0 || opt::mode == 1)
									break;
							}
						}
					}
				}

				bool skipped_repeat = false;
				string prev_insertion;
				// make our edit
				switch (best_edit_type) {
					case 1: // SUBSTITUTION MADE
						// apply the change to the actual contigSequence
						if (tNode->node_type == 0) {
							contigSeq[t_seq_i] = best_sub_base; 
							sRec subst; 
							subst.draft_char = draft_char;
							subst.pos = t_seq_i; 
							subst.sub_base = best_sub_base; 
							subst.num_support = best_num_support;
							substitution_record.push(subst); 
						} else if (tNode->node_type == 1)
							tNode->c = best_sub_base;
						// make sure we change our current hash to match it
						NTMC64_changelast(draft_char, best_sub_base, opt::k, opt::h, fhVal, rhVal, hVal); 
						if (opt::verbose) 
							std::cout << "\tt_seq_i: " << t_seq_i << " SUB: " << best_sub_base
								<< " check_present: " << best_num_support << std::endl;
						break;
					case 2: // INSERTION MADE
						// check if we need to check the insertion 
						// 	or low complexity and make that check before preceding
						prev_insertion = getPrevInsertion(t_seq_i, tNode);
						if (prev_insertion.size()+best_indel.size() >= opt::k) {
							std::cout << "\tprev_insertion: " << prev_insertion << std::endl;
							RopeLink *to_remove;
							if (isRepeatInsertion(prev_insertion)){
								std::cout << "\t\t is a repeat insertion" << std::endl;
								if (tNode != nullptr && tNode->node_type == 0 
										&& t_seq_i == tNode->s_pos)
								       	tNode = tNode->left;
								while (tNode != nullptr && tNode->node_type == 1) {
									to_remove = tNode; 
									std::cout << to_remove->node_type << " " << to_remove->c << std::endl; 
									tNode = tNode->left; 
									removeRopeLink(to_remove); 
								}
								NTMC64(findAcceptedKmer(h_seq_i, t_seq_i,
										       	hNode, tNode, contigSeq).c_str(), 
										opt::k, opt::h, fhVal, rhVal, hVal);
							        std::cout << "\tremoved prev_insertion: " << h_seq_i
									<< " " << t_seq_i << std::endl; 	
								skipped_repeat=true;
							}
							for (unsigned j=0; j<best_indel.size() && !skipped_repeat; j++) {
								prev_insertion += best_indel[j]; 
								std::cout << "\tnew_insertion: " << prev_insertion << std::endl; 
								if (isRepeatInsertion(prev_insertion)) {
									std::cout << "\t\t is a repeat insertion" << std::endl; 
									if (tNode != nullptr && tNode->node_type == 0
											&& t_seq_i == tNode->s_pos) 
										tNode = tNode->left; 	
									while (tNode != nullptr && tNode->node_type == 1) {
										to_remove = tNode; 
										tNode = tNode->left; 
										removeRopeLink(to_remove); 
									}
									NTMC64(findAcceptedKmer(h_seq_i, t_seq_i, 
												hNode, tNode, contigSeq).c_str(),
											opt::k, opt::h, fhVal, rhVal, hVal); 
									skipped_repeat=true;
								}	
							}
						}
						// if we didn't skip this region for a repeat, then make this insertion
						if (!skipped_repeat) {
							// make the insertion
							for (int j=best_indel.size()-1; j>=0; j--) {
								makeCharacterInsertion(tNode, t_seq_i, 
										best_indel[j], best_num_support);
							}
							NTMC64_changelast(draft_char, best_indel[0], 
									opt::k, opt::h, fhVal, rhVal, hVal); 
							if (opt::verbose) 
								std::cout << "\tt_seq_i: " << t_seq_i 
									<< " INS: " << best_indel << " check_present: " 
									<< best_num_support << std::endl; 
						}
						break;
					case 3: // DELETION MADE
						if (opt::verbose)
							std::cout << "\tt_seq_i: " << t_seq_i 
								<< " DEL: " << best_indel << " check_present: "	
								<< best_num_support << std::endl;
						makeDeletion(tNode, t_seq_i, best_indel.size(), best_num_support); 
						NTMC64_changelast(draft_char, after_t.at(best_indel.size()), 
								opt::k, opt::h, fhVal, rhVal, hVal); 
						break;
					case 0:
						if (opt::verbose)
							std::cout << "\tt_seq_i: " << t_seq_i << " FIX NOT FOUND" << std::endl; 
						break;
				}
			}
		}
		// roll and skip over non-ATGC containing kmers
		int target_t_seq_i = -1; 
		do {
			if (roll(h_seq_i, t_seq_i, hNode, tNode, fhVal, rhVal, hVal, contigSeq, charIn)) {
				if (!isAcceptedBase(toupper(charIn))) target_t_seq_i = t_seq_i + opt::k; 
			} else {
				continue_edit = false;
				break;
			}
		} while (target_t_seq_i >= 0 && t_seq_i != target_t_seq_i); 
	} while (continue_edit); 

	// clean allocated memory for hash
	delete [] hVal; 
	hVal = nullptr;

#pragma omp critical(write)
	{
		// write this to file	
		writeEditsToFile(dfout, rfout, contigHdr, contigSeq, master_root, substitution_record); 
	}

	// clean allocated memory for RopeLink
	cleanRopeLinks(master_root);
	master_root = nullptr;
	tNode = nullptr; 
	hNode = nullptr;
}	

/* Read the contigs from the file and polish each contig. */
void readAndCorrect(BloomFilter& bloom) {
	// read file handle
	gzFile dfp; 
	dfp = gzopen(opt::draft_filename.c_str(), "r"); 
	kseq_t * seq = kseq_init(dfp); 
//	bool stop = false; 
	int num_contigs=0; 


	// outfile handles
	string d_filename = opt::outfile_prefix+"_edited.fa"; 
	string r_filename = opt::outfile_prefix+"_changes.tsv"; 
	FILE* dfout = fopen(d_filename.c_str(), "w"); 
	FILE* rfout = fopen(r_filename.c_str(), "w");

	fprintf(rfout, 
		"ID\tbpPosition+1\tOriginalBase\tNewBase Support %d-mer (out of %d)\tAlternateNewBase\tAlt.Support %d-mers\n",
			opt::k, 
			(opt::k / 3),
			opt::k); 

#pragma omp parallel shared(seq,dfout,rfout)
	{
		string contigHdr, contigSeq, contigName;
		bool stop = false;

		while(1) {
#pragma omp critical(reading)
			{
				if (!stop && kseq_read(seq) >= 0) {
					contigHdr = seq->name.s; 
					if (seq->comment.l) contigName = contigHdr + " " + seq->comment.s; 
					else contigName = contigHdr; 
					contigSeq = seq->seq.s; 
				} else
					stop = true; 
			}
			if (stop) 
				break; 
			else {
				unsigned seq_len = contigSeq.length();
				if (opt::verbose) std::cout << contigName << std::endl; 
				if (seq_len >= opt::min_contig_len) {
					kmerizeAndCorrect(contigName, contigSeq, seq_len, bloom, dfout, rfout); 
				}
#pragma omp atomic
				num_contigs++; 
				if (num_contigs % 1000000 == 0) {
					std::cout << "Processed " << num_contigs << std::endl; 
				}
			}
		}
	}
//#pragma omp barrier
	kseq_destroy(seq); 
	gzclose(dfp); 
	fclose(dfout); 
	fclose(rfout); 
}

int main (int argc, char ** argv) {
	bool die = false; 
	for (int c; (c=getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
		std::istringstream arg(optarg != NULL ? optarg: ""); 
		switch (c) {
			case '?':
				die=true; 
				break; 
			case 't':
				arg >> opt::nthreads; 
				break;
			case 'f': 
				arg >> opt::draft_filename; 
				break; 
			case 'k': 
				arg >> opt::k; 
				break; 
			case 'z':
				arg >> opt::min_contig_len; 
				break; 
			case 'b':
				arg >> opt::outfile_prefix; 
				break;
			case 'r':
				arg >> opt::bloom_filename;
				break; 
			case 'd':
				arg >> opt::max_deletions;
				break; 
			case 'i':
				arg >> opt::max_insertions;
				break;
			case 'x':
				arg >> opt::missing_threshold;
				break;
			case 'y':
				arg >> opt::edit_threshold;
				break;
			case 'm':
				arg >> opt::mode; 
				break;
			case 'v':
				arg >> opt::verbose; 
				break;
			case OPT_HELP:
				std::cerr << USAGE_MESSAGE; 
				exit(EXIT_SUCCESS); 
			case OPT_VERSION:
				std::cerr << VERSION_MESSAGE; 
				exit(EXIT_SUCCESS); 
		}
		if (optarg != NULL && (!arg.eof() || arg.fail())) {
			std::cerr << PROGRAM ": invalid option: `-"
				<< (char)c << optarg << "'\n";
			exit(EXIT_FAILURE); 
		}
	}

	time_t rawtime; 
	time(&rawtime); 
	std::cout << "\n-----------------Running ntEdit------------- " << ctime(&rawtime); 

	// check the draft file is specified
	if (opt::draft_filename.empty()) {
		std::cerr << PROGRAM ": error: need to specify assembly draft file (-f)\n";
		die=true; 
	} else {
		// if the file is specified check that it is readable
		assert_readable(opt::draft_filename); 
	}

	// check that the bloom filter file is specified
	if (opt::bloom_filename.empty()) {
		std::cerr << PROGRAM ": error: need to specify the bloom filter file (-r)\n"; 
		die=true; 
	} else {
		// if the file is specified check that it is readable
		assert_readable(opt::bloom_filename); 
	}


	if (die) {
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n"; 
		exit(EXIT_FAILURE);
	}

	// check that the parameters x and y are in bound
	if (opt::missing_threshold < 3 && opt::missing_threshold > opt::k && 
			opt::edit_threshold < 3 && opt::edit_threshold > opt::k) {
		std::cerr << PROGRAM ": warning: x and y parameters must be >=3 and <=k; x and y were reset to default values x=5, y=9.\n"; 
		opt::missing_threshold = 5; 
		opt::missing_threshold = 5; 
		opt::edit_threshold=9; 
	}

	// check that the parameters i and d are in bound
	if ((opt::max_insertions == 0 && opt::max_deletions > 0) || (opt::max_insertions == 1 && opt::max_deletions > 1)) {
		std::cerr << PROGRAM ": warning: i and d parameter combination is not possible; d was set to the value of i.\n"; 
		opt::max_deletions = opt::max_insertions; 
	}


	// get the basename for the file
	string draft_basename = opt::draft_filename.substr(opt::draft_filename.find_last_of("/\\")+1); 
	string	bloom_basename = opt::bloom_filename.substr(opt::bloom_filename.find_last_of("/\\")+1); 

	// set the outfile prefix if it wasn't given 
	if (opt::outfile_prefix.empty()) {
		std::ostringstream outfile_name; 
		outfile_name << draft_basename
			<< "_k" << opt::k
			<< "_z" << opt::min_contig_len
			<< "_r" << bloom_basename
			<< "_i" << opt::max_insertions
			<< "_d" << opt::max_deletions
			<< "_m" << opt::mode;
		opt::outfile_prefix = outfile_name.str(); 
	}

	// print parameters: 
	std::cout << "Running: " << PROGRAM
		<< "\n -t " << opt::nthreads
		<< "\n -f " << draft_basename
		<< "\n -k " << opt::k
		<< "\n -z " << opt::min_contig_len
		<< "\n -b " << opt::outfile_prefix
		<< "\n -r " << bloom_basename
		<< "\n -i " << opt::max_insertions
		<< "\n -d " << opt::max_deletions
		<< "\n -m " << opt::mode
		<< "\n -v " << opt::verbose
		<< std::endl; 

	// Threading information 
	omp_set_num_threads(opt::nthreads); 

	// Load bloom filter
	time(&rawtime); 
	std::cout << "\n----------Loading Bloom Filter From File------ " << ctime(&rawtime);
	BloomFilter bloom(opt::bloom_filename.c_str());
	opt::h = bloom.getHashNum(); 

	// Checks for the bloom filter
	if (opt::h == 0) {
		std::cerr << PROGRAM ": error: Bloom Filter file is incorrect."; 
		exit(EXIT_FAILURE); 
	}

	if (opt::k != bloom.getKmerSize()) {
		std::cerr << PROGRAM ": error: Bloom Filter k size is different than ntEdit k size."; 
		exit(EXIT_FAILURE); 
	}	

	// print bloom filter details
	bloom.printBloomFilterDetails(); 
	
	// Read & edit contigs
	time(&rawtime); 
	std::cout << "\n----------Reading and Polishing Draft-------- " << ctime(&rawtime); 
	readAndCorrect(bloom); 

	time(&rawtime); 
	std::cout << "\n-----------ntEdit Polishing Complete!--------- " << ctime(&rawtime); 

	return 0; 
}
