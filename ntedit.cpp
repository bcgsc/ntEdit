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
#include <vector>
#include <unordered_set>
#include <unordered_map>
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
	PROGRAM " Version 0.0.1\n"
	"Written by Rene Warren, Hamid Mohamadi, and Jessica Zhang.\n"
	"Copyright 2018 Canada's Michael smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] = 
	"Usage: " PROGRAM " \n"
	"Make haploid edits to an assembly. \n"
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
	int verbose=0; 
}

static const char shortopts[] = "t:f:s:k:z:b:r:v:d:i:x:y:"; 

enum { OPT_HELP = 1, OPT_VERSION }; 

static const struct option longopts[] = {
	{"threads",	required_argument, NULL, 't'},
	{"draft_file",	required_argument, NULL, 'f'},
	{"k",	required_argument, NULL, 'k'},
	{"minimum contig length",	required_argument, NULL, 'z'},
	{"maximum insertions",	required_argument, NULL, 'i'},
	{"maximum deletions", required_argument, NULL, 'd'},
	{"edit threshold",	required_argument, NULL, 'y'},
	{"missing threshold",	required_argument, NULL, 'x'},
	{"bloom filename", required_argument, NULL, 'r'},
	{"outfile prefix", required_argument, NULL, 'b'},
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

struct FixInfo {
	unsigned pos; // 0-based index of fix; 
	unsigned char incorrect_base; 
	string indel; // always refers to what happens before or including t_seq_i
	unsigned char substitution; // always refers to the char at t_seq_i 
	unsigned num_support; 
};

// print an error message and exit if path is not readable
static inline void assert_readable(const string& path) {
	if (access(path.c_str(), R_OK) == -1) {
		std::cerr << PROGRAM ": error: `" << path << "': " << strerror(errno) << std::endl; 
		exit(EXIT_FAILURE); 
	}
}


unsigned countDelInRow(string indel, unsigned pos) {
	unsigned del_in_row=0;
	for (unsigned i=pos; i<indel.length(); i++) {
		if (indel[i]=='-') del_in_row++; 
		else break; ;
	}
	return del_in_row;
}

string getInsertionBlock(string indel, unsigned start_pos) {
	string insertion=""; 
	for (unsigned i=start_pos; i<indel.length(); i++) {
		if (isupper(indel[i])) insertion += indel[i]; 
		else break;
	}
	return insertion;
}

bool isAcceptedBase(unsigned char C) {
	return (C=='A' || C=='T' || C=='G' || C=='C'); 
}

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

bool getNextChars(unsigned char& charOut, unsigned char& charIn, unsigned& h_seq_i, unsigned& t_seq_i,
		int& num_indel_front, string& indel_bases,  
		int& num_indel_back, deque<unsigned char>& add_on_end,
		const string& contigSeq, unsigned seq_len) {

	// Deal with charIn
	//
	// If there is nothing to add on the back
	unsigned skipping_nonbase=0; 
	if (num_indel_back <= 0) {
		t_seq_i++; 
		if (t_seq_i>=seq_len) return false; 
		charIn=toupper(contigSeq.at(t_seq_i));
	        num_indel_back = 0; 	
		
	} else {
		if (num_indel_back==1) charIn=toupper(contigSeq.at(t_seq_i));
		else {
			charIn=add_on_end.front();
			add_on_end.pop_front();
		} 
		num_indel_back--;
	}

	if (num_indel_front==0) {
		h_seq_i++; 
	} else {
		charOut = indel_bases.at(indel_bases.length()-num_indel_front); 
		if (charOut == '-') {
			unsigned deletion_start = indel_bases.length()-num_indel_front; 
			unsigned num_del = countDelInRow(indel_bases, deletion_start); 
			h_seq_i += num_del; 
			//remove deletions so we don't do it again
			num_indel_front-=num_del;	
			charOut = toupper(contigSeq.at(h_seq_i));
		        h_seq_i++; 	
		} else num_indel_front--; 
	}

	return true; 
}

int getNumIndelFront(pair<FixInfo, bool>& edit) {
	if (!edit.second){
		edit.second=true; 
		return edit.first.indel.length(); 

	}	
	return 0;
}

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

bool isRepeatInsertion(string possible_repeat) {
	int n = possible_repeat.size(); 
	vector<int> lps(n); 

	computeLPSArray(possible_repeat, n, lps); 

	int len = lps[n-1]; 
	return (len>0 && n%(n-len) == 0);
}

void writeEditsToFile(FILE* dfout, FILE* rfout, vector<pair<FixInfo, bool>>& edits, 
			const string& contigHdr, const string& contigSeq) {
	fprintf(dfout, ">%s\n", contigHdr.c_str()); 
	for (unsigned i=0; i<contigSeq.length(); i++) {
		FixInfo fix_rec = edits[i].first; 
		bool del=false;
		for (unsigned j=0; j<fix_rec.indel.size(); j++) {
			if (isalpha(fix_rec.indel[j])) {
				string insertion=getInsertionBlock(fix_rec.indel, j);
				fprintf(dfout, "%s", insertion.c_str()); 
				fprintf(rfout, "%s\t%d\t%c\t%c%s\t%d\n", 
						contigHdr.c_str(),
						i+1,
						fix_rec.incorrect_base,
						'+',
						insertion.c_str(), 
						fix_rec.num_support);
				j += insertion.size()-1; 
			} else if (fix_rec.indel[j] == '-') {
				unsigned num_del = countDelInRow(fix_rec.indel, j); 
				string deletion;
				for (unsigned k=0; k<num_del; k++) {
					deletion += contigSeq[i+k]; 
				}
				fprintf(rfout, "%s\t%d\t%c\t%c%s\t%d\n", 
						contigHdr.c_str(),
						i+1,
						fix_rec.incorrect_base,
						'-',
						deletion.c_str(), 
						fix_rec.num_support);
				i+=num_del;
				j+=num_del; 
				del=true;
				break;
			}
		}
		if (!del && fix_rec.substitution) {
			fprintf(rfout, "%s\t%d\t%c\t%c\t%d\n", 
					contigHdr.c_str(),
					i+1,
					fix_rec.incorrect_base,
					fix_rec.substitution,
					fix_rec.num_support);
		}
		fprintf(dfout, "%c", contigSeq[i]);
	}
	fprintf(dfout, "\n");
}	

int tryDeletions(const unsigned char draft_char, unsigned& deletions_done, string& deleted_bases,
		unsigned& h_seq_i, unsigned& t_seq_i,
		uint64_t& fhVal, uint64_t& rhVal, uint64_t* hVal,
		int& num_indel_front, vector<pair<FixInfo, bool>>& edits, deque<unsigned>& subset_indel_pos,
		int& num_indel_back, deque<unsigned char>& add_on_end, 
		const string& contigSeq, unsigned seq_len, BloomFilter& bloom) {

	// update deletions done since we are doing one now
	deletions_done++; 

	// temporary values reset
	uint64_t temp_fhVal=fhVal;
	uint64_t temp_rhVal=rhVal;
	unsigned h_i=h_seq_i, t_i=t_seq_i; 
	int temp_num_indel_front=num_indel_front;
	int temp_num_indel_back=num_indel_back; 
	deque<unsigned char> temp_add_on_end=add_on_end; 
	unsigned char charOut, charIn;


	// make our deletion
	unsigned i=deletions_done; 
	unsigned char next_char;
	if (!temp_add_on_end.empty()) {
		next_char = temp_add_on_end.front(); 
		for (; i>0 && !temp_add_on_end.empty(); i--) {
			next_char = temp_add_on_end.front(); 
			temp_add_on_end.pop_front(); 
			temp_num_indel_back--; 
		}
	}
	if (i>0) {
		t_i += i-temp_num_indel_back;
	        temp_num_indel_back--; 	
		if (t_i>=seq_len) return 0; 
		next_char = contigSeq.at(t_i); 
		if (edits[t_i].first.indel.empty()) deleted_bases += toupper(contigSeq.at(t_i-1)); 
	}
	
	//track confirmations
	unsigned check_present=0;

	// update our hash
	NTMC64_changelast(draft_char, next_char, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
        if (bloom.contains(hVal)) check_present++; 	


	// check every 3rd kmer to confirm deletions
	for (unsigned k=1; k<=(opt::k-2) && h_i<seq_len; k++) {
		charOut=contigSeq.at(h_i); 
		if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
				&& !edits[h_i].second) { 
			temp_num_indel_front=getNumIndelFront(edits[h_i]);
			subset_indel_pos.push_back(h_i); 
		}
		if (getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front, edits[h_i].first.indel, 
					temp_num_indel_back, temp_add_on_end, contigSeq, seq_len)) {
			NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
			if (k%3==0 && bloom.contains(hVal)) check_present++; 
		}
	}

	// reset indels if we have changed it
	while (!subset_indel_pos.empty()) {
		edits[subset_indel_pos.front()].second=false; 
		subset_indel_pos.pop_front(); 
	}

	if (opt::verbose)
		std::cout << "\t\tdeleting: " << deleted_bases << " check_present: " << check_present << std::endl; 
	if (check_present>=((float)opt::k / opt::edit_threshold)) {
		unsigned leftover_del=deletions_done;
		if (!edits[t_seq_i].first.indel.empty() && num_indel_back > 0) {
			if (opt::verbose) std::cout << "\t\t\tprev edits[t_seq_i]: " << edits[t_seq_i].first.indel << std::endl; 
			if (deletions_done >= num_indel_back) {
				// deleted more than a previous insertion
				edits[t_seq_i].first.indel.erase(edits[t_seq_i].first.indel.length()-num_indel_back, num_indel_back); 
				// remove all of our insertions
				add_on_end.clear();
				leftover_del -= num_indel_back;
				// edge case where therew as a huge insertion before 
				if (h_seq_i == t_seq_i && num_indel_front > 0) {
					num_indel_front -= num_indel_back; 
				}
				num_indel_back=0; 
			} else {
				// did nto delete all of a previous insertion
				edits[t_seq_i].first.indel.erase(edits[t_seq_i].first.indel.length()-num_indel_back, deletions_done); 
				for (; leftover_del>0 && !add_on_end.empty(); leftover_del--) {
					add_on_end.pop_front();
					num_indel_back--;
				}
				// edge case where we had amde a insan insertion before
				if (h_seq_i==t_seq_i) {
					num_indel_front -= deletions_done; 
				}
			}
			if (opt::verbose) std::cout << "\t\t\tnew edits[t_seq_i]: " << edits[t_seq_i].first.indel << std::endl; 
		}
		if (h_seq_i == t_seq_i && num_indel_front>0) num_indel_front += leftover_del; 
		if (leftover_del > 0) {
			num_indel_back=(-1)*leftover_del; 
			string num_del_string=""; 
			for (unsigned j=0; j<leftover_del;j++) {
				num_del_string += '-'; 
			}
			edits[t_seq_i].first.incorrect_base = draft_char; 
			edits[t_seq_i].first.indel += num_del_string; 
			edits[t_seq_i].first.num_support = check_present; 
			num_indel_back =(-1)*leftover_del; 
			t_seq_i += leftover_del; 
		}

		// Update some things for next iteration 
		NTMC64_changelast(draft_char, next_char, opt::k, opt::h, fhVal, rhVal, hVal); 
		if (opt::verbose)
			std::cout << "\tt_seq_i: " << t_seq_i << " DEL: " << deleted_bases << " check_present "
				<< check_present << std::endl; 
		if (num_indel_back == 0) return -1; 
		return num_indel_back; 
	}

	return 0; 
}

int tryIndels(const unsigned char draft_char, const unsigned char index_char, unsigned& h_seq_i, unsigned& t_seq_i, 
		uint64_t& fhVal, uint64_t& rhVal, uint64_t* hVal, 
		int& num_indel_front, vector<pair<FixInfo, bool>>& edits, deque<unsigned>& subset_indel_pos,
		int& num_indel_back, deque<unsigned char>& add_on_end, 
		const string& contigSeq, unsigned seq_len, 
		BloomFilter& bloom, unsigned& deletions_done) {

	// for rec of deleted bases later
	string deleted_bases=""; 
	
	deque<unsigned char> temp_add_on_end;
	unsigned char charOut, charIn;
	
	// Try all of the combinations of indels starting with our index_char
	for (unsigned i=0; i<num_tries[opt::max_insertions]; i++) {

		// what are we inserting?
		string insertion_bases=multi_possible_bases[index_char][i]; 

		// temporary values reset
		uint64_t temp_fhVal=fhVal; 
		uint64_t temp_rhVal=rhVal; 
		unsigned h_i=h_seq_i, t_i=t_seq_i;
		int temp_num_indel_front=num_indel_front; 
		
		// temporary values specific for this insertion
		int temp_num_indel_back=insertion_bases.length(); 
		temp_add_on_end.clear(); 
		for(unsigned j=1; j<temp_num_indel_back; j++) {
			temp_add_on_end.push_back(insertion_bases[j]); 
		}

		// make first character insertion
		NTMC64_changelast(draft_char, index_char, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);

		// check every 3rd kmer to confirm insertion 
		unsigned check_present=0; 
		unsigned k=1;
		while (temp_num_indel_back>1 && h_i<seq_len) {
			charOut=contigSeq.at(h_i); 
			if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
					&& !edits[h_i].second) { 
				temp_num_indel_front=getNumIndelFront(edits[h_i]);
				subset_indel_pos.push_back(h_i); 
			}
			if (getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front, edits[h_i].first.indel, temp_num_indel_back, 
						temp_add_on_end, contigSeq, seq_len)) { 
				NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
				if (k%3==1 && bloom.contains(hVal)) check_present++; 
			}
			k++;
		}
		// do the original character
		if (h_i<seq_len) {
			charOut = contigSeq.at(h_i); 
			if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
					&& !edits[h_i].second) { 
				temp_num_indel_front=getNumIndelFront(edits[h_i]);
				subset_indel_pos.push_back(h_i); 
			}
			if (getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front, edits[h_i].first.indel, temp_num_indel_back,
						temp_add_on_end, contigSeq, seq_len)) {
				NTMC64(charOut, draft_char, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal); 
				if (k%3==1 && bloom.contains(hVal)) check_present++; 
			}
			k++; 
		}
		// temporary values specific for after insertion
		temp_num_indel_back=num_indel_back; 
		temp_add_on_end=add_on_end; 
		// check every 3rd kmer in rest of subset to confirm insertion
		for (; k<=(opt::k-1) && h_i<seq_len; k++) {
			charOut=contigSeq.at(h_i); 
			if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
					&& !edits[h_i].second) { 
				temp_num_indel_front=getNumIndelFront(edits[h_i]);
				subset_indel_pos.push_back(h_i); 
			}
			if (getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front, edits[h_i].first.indel, temp_num_indel_back,
						temp_add_on_end, contigSeq, seq_len)) {
				NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal); 
				if (k%3==1 && bloom.contains(hVal)) check_present++; 
			}
		}

		// reset indels if we have changed them
		while (!subset_indel_pos.empty()) {
			edits[subset_indel_pos.front()].second=false; 
			subset_indel_pos.pop_front(); 
		}

		if(opt::verbose) 
			std::cout << "\t\tinserting: " << insertion_bases << " check_present: " << check_present << std::endl;
		if (check_present>=((float) opt::k / opt::edit_threshold)) {
			if (edits[t_seq_i].first.indel.empty()) {
				//record for next iteration
				edits[t_seq_i].first.incorrect_base = draft_char; 
				edits[t_seq_i].first.indel = insertion_bases; 
				edits[t_seq_i].second=false; 
				num_indel_back=insertion_bases.length();
				add_on_end.clear(); 
				for (unsigned j=1; j<num_indel_back; j++) {
					add_on_end.push_back(insertion_bases[j]); 
				}
			} else {
				if (opt::verbose) std::cout << "\t\t\tprev insertion: " << edits[t_seq_i].first.indel;
				edits[t_seq_i].first.indel.insert(edits[t_seq_i].first.indel.size()-num_indel_back, insertion_bases); 
				if (edits[t_seq_i].first.indel.size() >= opt::k) {
					if (isRepeatInsertion(edits[t_seq_i].first.indel)) {
						// reset so that we can jump over it
						edits[t_seq_i].first.indel.clear();		
						insertion_bases = ""; 
						//num_indel_front=0; 
						num_indel_back = 0; 
						add_on_end.clear(); 
						// skip over this bad one and move on
						if (opt::verbose) std::cout 
							<< "\t\tremoved and jumped over low complexity this repeat insertion at: "
						       	<< t_seq_i << std::endl; 
						h_seq_i = findFirstAcceptedKmer(t_seq_i, contigSeq); 
						t_seq_i = h_seq_i+opt::k-1;
						if (h_seq_i < seq_len && t_seq_i < seq_len) {
							NTMC64(contigSeq.substr(h_seq_i, opt::k).c_str(), opt::k, opt::h, fhVal, rhVal, hVal); 
							charIn=contigSeq.at(t_seq_i); 
						}
						return num_indel_back; 
					}
				}
				if (opt::verbose) std::cout << " new insertion: " << edits[t_seq_i].first.indel << std::endl; 
				if (num_indel_back > 0) {
					add_on_end.push_front(draft_char); 	
					num_indel_back += insertion_bases.length();
				} else {
					num_indel_back += insertion_bases.length(); 
				}
				for (unsigned j=insertion_bases.length()-1; j>0; j--) {
					add_on_end.push_front(insertion_bases[j]); 
				}
			}
			// edge case for we weirdly made an insane insertion will need to adjust the front 
			if (h_seq_i==t_seq_i && num_indel_front > 0) {
				num_indel_front += insertion_bases.length();
			}
			edits[t_seq_i].first.num_support = check_present;

			// update some things for next iteration
			NTMC64_changelast(draft_char, index_char, opt::k, opt::h, fhVal, rhVal, hVal);
			if (opt::verbose)
				std::cout << "\tt_seq_i: " << t_seq_i << " INS: " << insertion_bases << " check_present: " 
					<< check_present << std::endl; 
			return num_indel_back;
		} else {
			// Try a deletion is there are any left
			if (deletions_done < opt::max_deletions) {
				int del_result= tryDeletions(draft_char, deletions_done, deleted_bases, h_seq_i, t_seq_i,
						fhVal, rhVal, hVal, num_indel_front, edits, 
						subset_indel_pos, num_indel_back,
						add_on_end, contigSeq, seq_len, bloom); 
				if (del_result!=0) return del_result;
			}
		}
	}
	return 0; 
}


void kmerizeAndCorrect(string& contigHdr, string& contigSeq, unsigned seq_len, BloomFilter& bloom, 
		FILE* dfout, FILE* rfout) {

	// initialize values for hashing
	uint64_t fhVal, rhVal;
	uint64_t* hVal = new uint64_t[opt::h]; 
	unsigned char charOut, charIn; 
	
	// initialize storages to track our changes
	vector<pair<FixInfo, bool>> edits(seq_len); // changed from unordered_map to vector for optimization purpsoes
	int num_indel_front=0, num_indel_back=0; 
	deque<unsigned char> add_on_end;

	// readjust the first character depending on the first N or nonATGC kmer
	unsigned h_seq_i = findFirstAcceptedKmer(0, contigSeq); 
	unsigned t_seq_i = h_seq_i+opt::k-1; 

	// intialize our seed kmer
	if (h_seq_i < seq_len && t_seq_i < seq_len) {
		NTMC64(contigSeq.substr(h_seq_i, opt::k).c_str(), opt::k, opt::h, fhVal, rhVal, hVal); 
		charIn=contigSeq.at(t_seq_i); 
	}

	while (h_seq_i < seq_len && t_seq_i < seq_len) {
		if (opt::verbose) 
			std::cout <<  h_seq_i << " " << t_seq_i << " " << charOut << " " << charIn << " " 
				<< hVal[0] << hVal[1] << hVal[2] << std::endl;
		
		if (!bloom.contains(hVal)) {

			// set temporary values
			uint64_t temp_fhVal=fhVal;
			uint64_t temp_rhVal=rhVal;
			unsigned h_i=h_seq_i; 
			unsigned t_i=t_seq_i; 
			int temp_num_indel_front=num_indel_front; 
			int temp_num_indel_back=num_indel_back; 
			deque<unsigned char> temp_add_on_end=add_on_end;

			// stack to keep track of changes to indels (the vector containing prev changes)
			deque<unsigned> subset_indel_pos;

			// draft_char
			unsigned char draft_char = toupper(charIn);	

			// Confirm that this kmer is missing by checking every 3rd kmer in subset
			unsigned check_missing=0;
			bool subset_contains_nonATGC=false;
			for (unsigned k=1; k<=opt::k && h_i<seq_len; k++) {
				charOut=contigSeq.at(h_i); 
				if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
						&& !edits[h_i].second) { 
					temp_num_indel_front=getNumIndelFront(edits[h_i]);
					subset_indel_pos.push_back(h_i); 
				}
				if (!getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front, edits[h_i].first.indel, 
							temp_num_indel_back, temp_add_on_end, contigSeq, seq_len)) break; 
			       	NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
				if (!isAcceptedBase(charIn)) {
					subset_contains_nonATGC=true;
					break;
				}
				if (k%3==1 && !bloom.contains(hVal)) {
					check_missing++; 
				}
			}

			// reset the indels that we might've changed
			while(!subset_indel_pos.empty()) {
				edits[subset_indel_pos.front()].second=false; 
				subset_indel_pos.pop_front(); 
			}


			if (opt::verbose)
				std::cout << "\tcheck_missing: " << check_missing << std::endl; 
			if (!subset_contains_nonATGC && check_missing>=((float) opt::k / opt::missing_threshold)) {
				// keeps track of what kind of changes we are trying to make
				bool substitutions_only=false; 
				unsigned deletions_done=0;
			        bool made_indel=false; 	

				// record the incorrect base
				edits[t_seq_i].first.pos = t_seq_i; 

				if (!isAcceptedBase(draft_char)) 
					draft_char='N';

				// SUBSTITUTIONS
				unsigned char best_sub_base; 
				unsigned best_sub_count=0; 
				for (const unsigned char sub_base : bases_array[draft_char]) {

					// temporary value reset
					temp_fhVal=fhVal;
					temp_rhVal=rhVal; 

					// make the subst change
					NTMC64_changelast(draft_char, sub_base, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal); 

					// only continue to pry if the substitution is found in bloom filter
					if (bloom.contains(hVal)) {

						// temporary value reset
						h_i=h_seq_i; 
						t_i=t_seq_i; 
						temp_num_indel_front=num_indel_front; 
						temp_num_indel_back=num_indel_back; 
						temp_add_on_end=add_on_end; 

						// check every 3rd kmer to confirm substitution
						unsigned check_present=0;
					        unsigned k=1;	
						for (; k<opt::k && h_i<seq_len; k++) {
							charOut=contigSeq.at(h_i); 
							if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
									&& !edits[h_i].second) { 
								temp_num_indel_front=getNumIndelFront(edits[h_i]);
								subset_indel_pos.push_back(h_i); 
							}
							if (getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front,
										edits[h_i].first.indel, temp_num_indel_back, 
										temp_add_on_end, contigSeq, seq_len)) {
								NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, 
										temp_rhVal, hVal); 
								if (k%3==1 && bloom.contains(hVal)) check_present++; 
							}
						}

						// do the last kmer
						if (k%3==1) {
							if (temp_num_indel_front==0 && edits[h_i].first.indel != ""
									&& !edits[h_i].second) { 
								temp_num_indel_front=getNumIndelFront(edits[h_i]);
								subset_indel_pos.push_back(h_i); 
							}
							if (getNextChars(charOut, charIn, h_i, t_i, temp_num_indel_front, 
										edits[h_i].first.indel, temp_num_indel_back, 
										temp_add_on_end, 
										contigSeq, seq_len)) {
								NTMC64(sub_base, charIn, opt::k, opt::h, temp_fhVal,
										temp_rhVal, hVal); 
								if (bloom.contains(hVal)) check_present++; 
							}
						}

						// reset the indels we mightve changed
						while (!subset_indel_pos.empty()) {
							edits[subset_indel_pos.front()].second=false; 
							subset_indel_pos.pop_front(); 
						}
						

						if (opt::verbose) 
							std::cout << "\t\tsub: " << sub_base << " check_present: " 
								<< check_present << std::endl; 
						// if it is good, only do substitutions
						if (check_present>=((float) opt::k / opt::edit_threshold)) {
							substitutions_only=true; 
							if (check_present>best_sub_count) {
								best_sub_base = sub_base; 
								best_sub_count=check_present; 
							}
						} else if (!substitutions_only) {
							if(tryIndels(draft_char, sub_base, h_seq_i, t_seq_i,
									fhVal, rhVal, hVal, num_indel_front, edits,
									subset_indel_pos,
									num_indel_back, add_on_end, contigSeq, seq_len, bloom, 
									deletions_done)!=0) { 
							       made_indel=true;
						       	       break; 
							}
						}
					}
				}

				if (substitutions_only) {
					// This means we found a good substitution
					if (edits[t_seq_i].first.indel.empty() || num_indel_back==0) {
						contigSeq[t_seq_i]=best_sub_base;
						// Record our change
						edits[t_seq_i].first.substitution = best_sub_base;
					        edits[t_seq_i].first.incorrect_base = draft_char; 	
					} else {
						//update it for later when we remove
						edits[t_seq_i].first.indel[edits[t_seq_i].first.indel.size()-num_indel_back] = best_sub_base;
					}
					edits[t_seq_i].first.num_support = best_sub_count; 
					NTMC64_changelast(draft_char, best_sub_base, opt::k, opt::h, fhVal, rhVal, hVal); 
					if (opt::verbose)
						std::cout << "\tt_seq_i: " << t_seq_i << " SUB: " << best_sub_base 
							<< " check_present: " << best_sub_count << std::endl; 
				} else if (!made_indel) {
					// There was no good change
					if (opt::verbose)
						std::cout << "\tt_seq_i: " << t_seq_i << " FIX NOT FOUND." << std::endl; 
				}
			}
		}
		// Update to the next kmer
		int target_t_seq_i=-1;
		do {
			charOut=toupper(contigSeq.at(h_seq_i));
			if (num_indel_front==0 && edits[h_seq_i].first.indel != "" && !edits[h_seq_i].second) {
				num_indel_front=getNumIndelFront(edits[h_seq_i]);
			}
			if (!getNextChars(charOut, charIn, h_seq_i, t_seq_i, num_indel_front, edits[h_seq_i].first.indel, num_indel_back, add_on_end, contigSeq, seq_len)) break; 
			NTMC64(charOut, charIn, opt::k, opt::h, fhVal, rhVal, hVal); 
			if (!isAcceptedBase(charIn)) target_t_seq_i = t_seq_i + opt::k;
		} while (target_t_seq_i >=0 && t_seq_i != target_t_seq_i); 
	}

#pragma omp critical(writing)
	{
		// write this to file	
		writeEditsToFile(dfout, rfout, edits, contigHdr, contigSeq); 
	}
}	

void readAndCorrect(BloomFilter& bloom) {
	// read file handle
	int l;
	gzFile dfp; 
	dfp = gzopen(opt::draft_filename.c_str(), "r"); 
	kseq_t * seq = kseq_init(dfp); 


	// outfile handles
	string d_filename = opt::outfile_prefix+"_edited.fa"; 
	string r_filename = opt::outfile_prefix+"_changes.tsv"; 
	FILE* dfout = fopen(d_filename.c_str(), "w"); 
	FILE* rfout = fopen(r_filename.c_str(), "w");

	fprintf(rfout, 
		"ID\tbpPosition+1\tOriginalBase\tNewBase Support %d-mer (out of %d)\tAlternateNewBase\tAlt.Support %d-mers\n",
			opt::k, 
			(unsigned) opt::edit_threshold,
			opt::k); 

#pragma omp parallel
	{
		string contigHdr, contigSeq, contigName;
		int num_contigs=0; 

		while(1) {
#pragma omp critical(reading)
			{
				l = kseq_read(seq); 
				contigHdr = seq->name.s;
				if (seq->comment.l) contigName = contigHdr + " " + seq->comment.s; 
				else contigName = contigHdr; 
				contigSeq = seq->seq.s; 
			}
			if (l<0) break; 

			unsigned seq_len = contigSeq.length();
			if (opt::verbose) 
				std::cout << contigName << std::endl; 
			if (seq_len >= opt::min_contig_len) {
				kmerizeAndCorrect(contigName, contigSeq, seq_len, bloom, dfout, rfout); 
			}
			num_contigs++; 
			if (num_contigs % 1000000 == 0) {
				std::cout << "Processed " << num_contigs << std::endl; 
			}
		}
	}
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
			<< "_d" << opt::max_deletions;
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
	std::cout << "\n----------Reading and Correcting Draft-------- " << ctime(&rawtime); 
	readAndCorrect(bloom); 

	time(&rawtime); 
	std::cout << "\n-----------ntEdit Polishing Complete!--------- " << ctime(&rawtime); 

	return 0; 
}
