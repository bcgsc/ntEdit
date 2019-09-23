#define PROGRAM "ntEdit" // NOLINT

// clang-format off
#include <iostream> //NOLINT(llvm-include-order)
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
#include "lib/nthash.hpp" // NOLINT
#include "lib/BloomFilter.hpp"
// clang-format on

// NOLINTNEXTLINE
KSEQ_INIT(gzFile, gzread)

// NOLINTNEXTLINE(modernize-avoid-c-arrays)
static const char VERSION_MESSAGE[] =
    PROGRAM " Version 1.2.3\n"
            "Written by Rene Warren, Hamid Mohamadi, and Jessica Zhang.\n"
            "Copyright 2018, 2019 Canada's Michael smith Genome Science Centre\n";

// NOLINTNEXTLINE(modernize-avoid-c-arrays)
static const char USAGE_MESSAGE[] = PROGRAM
    " v1.2.3\n"
    "\n"
    "Scalable genome sequence polishing.\n"
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
    "	-X, 	ratio of number of kmers in the k subset that should be missing in order to attempt fix (higher=stringent), [default=0.5]\n"
    "	-Y, 	ratio of number of kmers in the k subset that should be present to accept an edit (higher=stringent), [default=0.5]\n"
    "	-c,	cap for the number of base insertions that can be made at one position, [default=k*1.5]\n"
    "	-j, 	controls size of kmer subset. When checking subset of kmers, check every jth kmer, [default=3]\n"
    "	-m,	mode of editing, range 0-2, [default=0]\n"
    "			0: best substitution, or first good indel\n"
    "			1: best substitution, or best indel\n"
    "			2: best edit overall (suggestion that you reduce i and d for performance)\n"
    "	-v,	verbose mode (-v 1 = yes, default = 0, no)\n"
    "\n"
    "	--help,		display this message and exit \n"
    "	--version,	output version information and exit\n"
    "\n"
    "	If one of X/Y is set, ntEdit will use those parameters instead. Otherwise, it uses x/y by default.\n"
    "\n";


namespace opt {
/* Defining magical numbers. */
constexpr int default_min_contig_len = 100;
constexpr int default_max_insertions = 4;
constexpr int default_max_deletions = 5;
constexpr float default_edit_threshold = 9.0000;
constexpr float default_missing_threshold = 5.0000;
constexpr float default_insertion_cap_ratio = 1.5;
float edit_ratio = 0.5;
float missing_ratio = 0.5;
bool use_ratio = false;
unsigned jump = 3;
unsigned nthreads = 1;
std::string draft_filename; // NOLINT
std::string bloom_filename; // NOLINT
std::string outfile_prefix; // NOLINT
unsigned k;
unsigned h = 0;
unsigned min_contig_len = default_min_contig_len;
unsigned max_insertions = default_max_insertions;
unsigned max_deletions = default_max_deletions;
float edit_threshold = default_edit_threshold;
float missing_threshold = default_missing_threshold;
unsigned insertion_cap = static_cast<unsigned>(static_cast<float>(opt::k) * default_insertion_cap_ratio);
int mode = 0;
int verbose = 0;
} // namespace opt

static const char shortopts[] = "t:f:s:k:z:b:r:v:d:i:X:Y:x:y:m:c:j:";

enum
{
	OPT_HELP = 1,
	OPT_VERSION
};

static const struct option longopts[] = {
	{ "threads", required_argument, nullptr, 't' },
	{ "draft_file", required_argument, nullptr, 'f' },
	{ "k", required_argument, nullptr, 'k' },
	{ "minimum_contig_length", required_argument, nullptr, 'z' },
	{ "maximum_insertions", required_argument, nullptr, 'i' },
	{ "maximum_deletions", required_argument, nullptr, 'd' },
	{ "insertion_cap", required_argument, nullptr, 'c' },
	{ "edit_threshold", required_argument, nullptr, 'y' },
	{ "missing_threshold", required_argument, nullptr, 'x' },
	{ "edit_ratio", required_argument, NULL, 'Y' },
	{ "missing_ratio", required_argument, NULL, 'X' },
	{ "jump", required_argument, NULL, 'j' },
	{ "bloom_filename", required_argument, nullptr, 'r' },
	{ "outfile_prefix", required_argument, nullptr, 'b' },
	{ "mode", required_argument, nullptr, 'm' },
	{ "verbose", required_argument, nullptr, 'v' },
	{ "help", no_argument, nullptr, OPT_HELP },
	{ "version", no_argument, nullptr, OPT_VERSION },
	{ nullptr, 0, nullptr, 0 }
};

// Setting up the number of tries when for each number of base insertion
std::vector<int> num_tries = { 0, 1, 5, 21, 85, 341 }; // NOLINT

// Setting up base array
// NOLINTNEXTLINE
std::unordered_map<unsigned char, std::vector<unsigned char>> bases_array = {
	{ 'A', { 'T', 'C', 'G' } },
	{ 'T', { 'A', 'C', 'G' } },
	{ 'C', { 'A', 'T', 'G' } },
	{ 'G', { 'A', 'T', 'C' } },
	{ 'N', { 'A', 'T', 'C', 'G' } }
};

// Setting all the indel combos
// NOLINTNEXTLINE
std::unordered_map<unsigned char, std::vector<string>> multi_possible_bases = {
	{ 'A',
	  { "A",     "AA",    "AC",    "AG",    "AT",    "AAA",   "AAC",   "AAG",   "AAT",   "ACA",
	    "ACC",   "ACG",   "ACT",   "AGA",   "AGC",   "AGG",   "AGT",   "ATA",   "ATC",   "ATG",
	    "ATT",   "AAAA",  "AAAC",  "AAAG",  "AAAT",  "AACA",  "AACC",  "AACG",  "AACT",  "AAGA",
	    "AAGC",  "AAGG",  "AAGT",  "AATA",  "AATC",  "AATG",  "AATT",  "ACAA",  "ACAC",  "ACAG",
	    "ACAT",  "ACCA",  "ACCC",  "ACCG",  "ACCT",  "ACGA",  "ACGC",  "ACGG",  "ACGT",  "ACTA",
	    "ACTC",  "ACTG",  "ACTT",  "AGAA",  "AGAC",  "AGAG",  "AGAT",  "AGCA",  "AGCC",  "AGCG",
	    "AGCT",  "AGGA",  "AGGC",  "AGGG",  "AGGT",  "AGTA",  "AGTC",  "AGTG",  "AGTT",  "ATAA",
	    "ATAC",  "ATAG",  "ATAT",  "ATCA",  "ATCC",  "ATCG",  "ATCT",  "ATGA",  "ATGC",  "ATGG",
	    "ATGT",  "ATTA",  "ATTC",  "ATTG",  "ATTT",  "AAAAA", "AAAAC", "AAAAG", "AAAAT", "AAACA",
	    "AAACC", "AAACG", "AAACT", "AAAGA", "AAAGC", "AAAGG", "AAAGT", "AAATA", "AAATC", "AAATG",
	    "AAATT", "AACAA", "AACAC", "AACAG", "AACAT", "AACCA", "AACCC", "AACCG", "AACCT", "AACGA",
	    "AACGC", "AACGG", "AACGT", "AACTA", "AACTC", "AACTG", "AACTT", "AAGAA", "AAGAC", "AAGAG",
	    "AAGAT", "AAGCA", "AAGCC", "AAGCG", "AAGCT", "AAGGA", "AAGGC", "AAGGG", "AAGGT", "AAGTA",
	    "AAGTC", "AAGTG", "AAGTT", "AATAA", "AATAC", "AATAG", "AATAT", "AATCA", "AATCC", "AATCG",
	    "AATCT", "AATGA", "AATGC", "AATGG", "AATGT", "AATTA", "AATTC", "AATTG", "AATTT", "ACAAA",
	    "ACAAC", "ACAAG", "ACAAT", "ACACA", "ACACC", "ACACG", "ACACT", "ACAGA", "ACAGC", "ACAGG",
	    "ACAGT", "ACATA", "ACATC", "ACATG", "ACATT", "ACCAA", "ACCAC", "ACCAG", "ACCAT", "ACCCA",
	    "ACCCC", "ACCCG", "ACCCT", "ACCGA", "ACCGC", "ACCGG", "ACCGT", "ACCTA", "ACCTC", "ACCTG",
	    "ACCTT", "ACGAA", "ACGAC", "ACGAG", "ACGAT", "ACGCA", "ACGCC", "ACGCG", "ACGCT", "ACGGA",
	    "ACGGC", "ACGGG", "ACGGT", "ACGTA", "ACGTC", "ACGTG", "ACGTT", "ACTAA", "ACTAC", "ACTAG",
	    "ACTAT", "ACTCA", "ACTCC", "ACTCG", "ACTCT", "ACTGA", "ACTGC", "ACTGG", "ACTGT", "ACTTA",
	    "ACTTC", "ACTTG", "ACTTT", "AGAAA", "AGAAC", "AGAAG", "AGAAT", "AGACA", "AGACC", "AGACG",
	    "AGACT", "AGAGA", "AGAGC", "AGAGG", "AGAGT", "AGATA", "AGATC", "AGATG", "AGATT", "AGCAA",
	    "AGCAC", "AGCAG", "AGCAT", "AGCCA", "AGCCC", "AGCCG", "AGCCT", "AGCGA", "AGCGC", "AGCGG",
	    "AGCGT", "AGCTA", "AGCTC", "AGCTG", "AGCTT", "AGGAA", "AGGAC", "AGGAG", "AGGAT", "AGGCA",
	    "AGGCC", "AGGCG", "AGGCT", "AGGGA", "AGGGC", "AGGGG", "AGGGT", "AGGTA", "AGGTC", "AGGTG",
	    "AGGTT", "AGTAA", "AGTAC", "AGTAG", "AGTAT", "AGTCA", "AGTCC", "AGTCG", "AGTCT", "AGTGA",
	    "AGTGC", "AGTGG", "AGTGT", "AGTTA", "AGTTC", "AGTTG", "AGTTT", "ATAAA", "ATAAC", "ATAAG",
	    "ATAAT", "ATACA", "ATACC", "ATACG", "ATACT", "ATAGA", "ATAGC", "ATAGG", "ATAGT", "ATATA",
	    "ATATC", "ATATG", "ATATT", "ATCAA", "ATCAC", "ATCAG", "ATCAT", "ATCCA", "ATCCC", "ATCCG",
	    "ATCCT", "ATCGA", "ATCGC", "ATCGG", "ATCGT", "ATCTA", "ATCTC", "ATCTG", "ATCTT", "ATGAA",
	    "ATGAC", "ATGAG", "ATGAT", "ATGCA", "ATGCC", "ATGCG", "ATGCT", "ATGGA", "ATGGC", "ATGGG",
	    "ATGGT", "ATGTA", "ATGTC", "ATGTG", "ATGTT", "ATTAA", "ATTAC", "ATTAG", "ATTAT", "ATTCA",
	    "ATTCC", "ATTCG", "ATTCT", "ATTGA", "ATTGC", "ATTGG", "ATTGT", "ATTTA", "ATTTC", "ATTTG",
	    "ATTTT" } },
	{ 'C',
	  { "C",     "CA",    "CC",    "CG",    "CT",    "CAA",   "CAC",   "CAG",   "CAT",   "CCA",
	    "CCC",   "CCG",   "CCT",   "CGA",   "CGC",   "CGG",   "CGT",   "CTA",   "CTC",   "CTG",
	    "CTT",   "CAAA",  "CAAC",  "CAAG",  "CAAT",  "CACA",  "CACC",  "CACG",  "CACT",  "CAGA",
	    "CAGC",  "CAGG",  "CAGT",  "CATA",  "CATC",  "CATG",  "CATT",  "CCAA",  "CCAC",  "CCAG",
	    "CCAT",  "CCCA",  "CCCC",  "CCCG",  "CCCT",  "CCGA",  "CCGC",  "CCGG",  "CCGT",  "CCTA",
	    "CCTC",  "CCTG",  "CCTT",  "CGAA",  "CGAC",  "CGAG",  "CGAT",  "CGCA",  "CGCC",  "CGCG",
	    "CGCT",  "CGGA",  "CGGC",  "CGGG",  "CGGT",  "CGTA",  "CGTC",  "CGTG",  "CGTT",  "CTAA",
	    "CTAC",  "CTAG",  "CTAT",  "CTCA",  "CTCC",  "CTCG",  "CTCT",  "CTGA",  "CTGC",  "CTGG",
	    "CTGT",  "CTTA",  "CTTC",  "CTTG",  "CTTT",  "CAAAA", "CAAAC", "CAAAG", "CAAAT", "CAACA",
	    "CAACC", "CAACG", "CAACT", "CAAGA", "CAAGC", "CAAGG", "CAAGT", "CAATA", "CAATC", "CAATG",
	    "CAATT", "CACAA", "CACAC", "CACAG", "CACAT", "CACCA", "CACCC", "CACCG", "CACCT", "CACGA",
	    "CACGC", "CACGG", "CACGT", "CACTA", "CACTC", "CACTG", "CACTT", "CAGAA", "CAGAC", "CAGAG",
	    "CAGAT", "CAGCA", "CAGCC", "CAGCG", "CAGCT", "CAGGA", "CAGGC", "CAGGG", "CAGGT", "CAGTA",
	    "CAGTC", "CAGTG", "CAGTT", "CATAA", "CATAC", "CATAG", "CATAT", "CATCA", "CATCC", "CATCG",
	    "CATCT", "CATGA", "CATGC", "CATGG", "CATGT", "CATTA", "CATTC", "CATTG", "CATTT", "CCAAA",
	    "CCAAC", "CCAAG", "CCAAT", "CCACA", "CCACC", "CCACG", "CCACT", "CCAGA", "CCAGC", "CCAGG",
	    "CCAGT", "CCATA", "CCATC", "CCATG", "CCATT", "CCCAA", "CCCAC", "CCCAG", "CCCAT", "CCCCA",
	    "CCCCC", "CCCCG", "CCCCT", "CCCGA", "CCCGC", "CCCGG", "CCCGT", "CCCTA", "CCCTC", "CCCTG",
	    "CCCTT", "CCGAA", "CCGAC", "CCGAG", "CCGAT", "CCGCA", "CCGCC", "CCGCG", "CCGCT", "CCGGA",
	    "CCGGC", "CCGGG", "CCGGT", "CCGTA", "CCGTC", "CCGTG", "CCGTT", "CCTAA", "CCTAC", "CCTAG",
	    "CCTAT", "CCTCA", "CCTCC", "CCTCG", "CCTCT", "CCTGA", "CCTGC", "CCTGG", "CCTGT", "CCTTA",
	    "CCTTC", "CCTTG", "CCTTT", "CGAAA", "CGAAC", "CGAAG", "CGAAT", "CGACA", "CGACC", "CGACG",
	    "CGACT", "CGAGA", "CGAGC", "CGAGG", "CGAGT", "CGATA", "CGATC", "CGATG", "CGATT", "CGCAA",
	    "CGCAC", "CGCAG", "CGCAT", "CGCCA", "CGCCC", "CGCCG", "CGCCT", "CGCGA", "CGCGC", "CGCGG",
	    "CGCGT", "CGCTA", "CGCTC", "CGCTG", "CGCTT", "CGGAA", "CGGAC", "CGGAG", "CGGAT", "CGGCA",
	    "CGGCC", "CGGCG", "CGGCT", "CGGGA", "CGGGC", "CGGGG", "CGGGT", "CGGTA", "CGGTC", "CGGTG",
	    "CGGTT", "CGTAA", "CGTAC", "CGTAG", "CGTAT", "CGTCA", "CGTCC", "CGTCG", "CGTCT", "CGTGA",
	    "CGTGC", "CGTGG", "CGTGT", "CGTTA", "CGTTC", "CGTTG", "CGTTT", "CTAAA", "CTAAC", "CTAAG",
	    "CTAAT", "CTACA", "CTACC", "CTACG", "CTACT", "CTAGA", "CTAGC", "CTAGG", "CTAGT", "CTATA",
	    "CTATC", "CTATG", "CTATT", "CTCAA", "CTCAC", "CTCAG", "CTCAT", "CTCCA", "CTCCC", "CTCCG",
	    "CTCCT", "CTCGA", "CTCGC", "CTCGG", "CTCGT", "CTCTA", "CTCTC", "CTCTG", "CTCTT", "CTGAA",
	    "CTGAC", "CTGAG", "CTGAT", "CTGCA", "CTGCC", "CTGCG", "CTGCT", "CTGGA", "CTGGC", "CTGGG",
	    "CTGGT", "CTGTA", "CTGTC", "CTGTG", "CTGTT", "CTTAA", "CTTAC", "CTTAG", "CTTAT", "CTTCA",
	    "CTTCC", "CTTCG", "CTTCT", "CTTGA", "CTTGC", "CTTGG", "CTTGT", "CTTTA", "CTTTC", "CTTTG",
	    "CTTTT" } },
	{ 'G',
	  { "G",     "GA",    "GC",    "GG",    "GT",    "GAA",   "GAC",   "GAG",   "GAT",   "GCA",
	    "GCC",   "GCG",   "GCT",   "GGA",   "GGC",   "GGG",   "GGT",   "GTA",   "GTC",   "GTG",
	    "GTT",   "GAAA",  "GAAC",  "GAAG",  "GAAT",  "GACA",  "GACC",  "GACG",  "GACT",  "GAGA",
	    "GAGC",  "GAGG",  "GAGT",  "GATA",  "GATC",  "GATG",  "GATT",  "GCAA",  "GCAC",  "GCAG",
	    "GCAT",  "GCCA",  "GCCC",  "GCCG",  "GCCT",  "GCGA",  "GCGC",  "GCGG",  "GCGT",  "GCTA",
	    "GCTC",  "GCTG",  "GCTT",  "GGAA",  "GGAC",  "GGAG",  "GGAT",  "GGCA",  "GGCC",  "GGCG",
	    "GGCT",  "GGGA",  "GGGC",  "GGGG",  "GGGT",  "GGTA",  "GGTC",  "GGTG",  "GGTT",  "GTAA",
	    "GTAC",  "GTAG",  "GTAT",  "GTCA",  "GTCC",  "GTCG",  "GTCT",  "GTGA",  "GTGC",  "GTGG",
	    "GTGT",  "GTTA",  "GTTC",  "GTTG",  "GTTT",  "GAAAA", "GAAAC", "GAAAG", "GAAAT", "GAACA",
	    "GAACC", "GAACG", "GAACT", "GAAGA", "GAAGC", "GAAGG", "GAAGT", "GAATA", "GAATC", "GAATG",
	    "GAATT", "GACAA", "GACAC", "GACAG", "GACAT", "GACCA", "GACCC", "GACCG", "GACCT", "GACGA",
	    "GACGC", "GACGG", "GACGT", "GACTA", "GACTC", "GACTG", "GACTT", "GAGAA", "GAGAC", "GAGAG",
	    "GAGAT", "GAGCA", "GAGCC", "GAGCG", "GAGCT", "GAGGA", "GAGGC", "GAGGG", "GAGGT", "GAGTA",
	    "GAGTC", "GAGTG", "GAGTT", "GATAA", "GATAC", "GATAG", "GATAT", "GATCA", "GATCC", "GATCG",
	    "GATCT", "GATGA", "GATGC", "GATGG", "GATGT", "GATTA", "GATTC", "GATTG", "GATTT", "GCAAA",
	    "GCAAC", "GCAAG", "GCAAT", "GCACA", "GCACC", "GCACG", "GCACT", "GCAGA", "GCAGC", "GCAGG",
	    "GCAGT", "GCATA", "GCATC", "GCATG", "GCATT", "GCCAA", "GCCAC", "GCCAG", "GCCAT", "GCCCA",
	    "GCCCC", "GCCCG", "GCCCT", "GCCGA", "GCCGC", "GCCGG", "GCCGT", "GCCTA", "GCCTC", "GCCTG",
	    "GCCTT", "GCGAA", "GCGAC", "GCGAG", "GCGAT", "GCGCA", "GCGCC", "GCGCG", "GCGCT", "GCGGA",
	    "GCGGC", "GCGGG", "GCGGT", "GCGTA", "GCGTC", "GCGTG", "GCGTT", "GCTAA", "GCTAC", "GCTAG",
	    "GCTAT", "GCTCA", "GCTCC", "GCTCG", "GCTCT", "GCTGA", "GCTGC", "GCTGG", "GCTGT", "GCTTA",
	    "GCTTC", "GCTTG", "GCTTT", "GGAAA", "GGAAC", "GGAAG", "GGAAT", "GGACA", "GGACC", "GGACG",
	    "GGACT", "GGAGA", "GGAGC", "GGAGG", "GGAGT", "GGATA", "GGATC", "GGATG", "GGATT", "GGCAA",
	    "GGCAC", "GGCAG", "GGCAT", "GGCCA", "GGCCC", "GGCCG", "GGCCT", "GGCGA", "GGCGC", "GGCGG",
	    "GGCGT", "GGCTA", "GGCTC", "GGCTG", "GGCTT", "GGGAA", "GGGAC", "GGGAG", "GGGAT", "GGGCA",
	    "GGGCC", "GGGCG", "GGGCT", "GGGGA", "GGGGC", "GGGGG", "GGGGT", "GGGTA", "GGGTC", "GGGTG",
	    "GGGTT", "GGTAA", "GGTAC", "GGTAG", "GGTAT", "GGTCA", "GGTCC", "GGTCG", "GGTCT", "GGTGA",
	    "GGTGC", "GGTGG", "GGTGT", "GGTTA", "GGTTC", "GGTTG", "GGTTT", "GTAAA", "GTAAC", "GTAAG",
	    "GTAAT", "GTACA", "GTACC", "GTACG", "GTACT", "GTAGA", "GTAGC", "GTAGG", "GTAGT", "GTATA",
	    "GTATC", "GTATG", "GTATT", "GTCAA", "GTCAC", "GTCAG", "GTCAT", "GTCCA", "GTCCC", "GTCCG",
	    "GTCCT", "GTCGA", "GTCGC", "GTCGG", "GTCGT", "GTCTA", "GTCTC", "GTCTG", "GTCTT", "GTGAA",
	    "GTGAC", "GTGAG", "GTGAT", "GTGCA", "GTGCC", "GTGCG", "GTGCT", "GTGGA", "GTGGC", "GTGGG",
	    "GTGGT", "GTGTA", "GTGTC", "GTGTG", "GTGTT", "GTTAA", "GTTAC", "GTTAG", "GTTAT", "GTTCA",
	    "GTTCC", "GTTCG", "GTTCT", "GTTGA", "GTTGC", "GTTGG", "GTTGT", "GTTTA", "GTTTC", "GTTTG",
	    "GTTTT" } },
	{ 'T',
	  { "T",     "TA",    "TC",    "TG",    "TT",    "TAA",   "TAC",   "TAG",   "TAT",   "TCA",
	    "TCC",   "TCG",   "TCT",   "TGA",   "TGC",   "TGG",   "TGT",   "TTA",   "TTC",   "TTG",
	    "TTT",   "TAAA",  "TAAC",  "TAAG",  "TAAT",  "TACA",  "TACC",  "TACG",  "TACT",  "TAGA",
	    "TAGC",  "TAGG",  "TAGT",  "TATA",  "TATC",  "TATG",  "TATT",  "TCAA",  "TCAC",  "TCAG",
	    "TCAT",  "TCCA",  "TCCC",  "TCCG",  "TCCT",  "TCGA",  "TCGC",  "TCGG",  "TCGT",  "TCTA",
	    "TCTC",  "TCTG",  "TCTT",  "TGAA",  "TGAC",  "TGAG",  "TGAT",  "TGCA",  "TGCC",  "TGCG",
	    "TGCT",  "TGGA",  "TGGC",  "TGGG",  "TGGT",  "TGTA",  "TGTC",  "TGTG",  "TGTT",  "TTAA",
	    "TTAC",  "TTAG",  "TTAT",  "TTCA",  "TTCC",  "TTCG",  "TTCT",  "TTGA",  "TTGC",  "TTGG",
	    "TTGT",  "TTTA",  "TTTC",  "TTTG",  "TTTT",  "TAAAA", "TAAAC", "TAAAG", "TAAAT", "TAACA",
	    "TAACC", "TAACG", "TAACT", "TAAGA", "TAAGC", "TAAGG", "TAAGT", "TAATA", "TAATC", "TAATG",
	    "TAATT", "TACAA", "TACAC", "TACAG", "TACAT", "TACCA", "TACCC", "TACCG", "TACCT", "TACGA",
	    "TACGC", "TACGG", "TACGT", "TACTA", "TACTC", "TACTG", "TACTT", "TAGAA", "TAGAC", "TAGAG",
	    "TAGAT", "TAGCA", "TAGCC", "TAGCG", "TAGCT", "TAGGA", "TAGGC", "TAGGG", "TAGGT", "TAGTA",
	    "TAGTC", "TAGTG", "TAGTT", "TATAA", "TATAC", "TATAG", "TATAT", "TATCA", "TATCC", "TATCG",
	    "TATCT", "TATGA", "TATGC", "TATGG", "TATGT", "TATTA", "TATTC", "TATTG", "TATTT", "TCAAA",
	    "TCAAC", "TCAAG", "TCAAT", "TCACA", "TCACC", "TCACG", "TCACT", "TCAGA", "TCAGC", "TCAGG",
	    "TCAGT", "TCATA", "TCATC", "TCATG", "TCATT", "TCCAA", "TCCAC", "TCCAG", "TCCAT", "TCCCA",
	    "TCCCC", "TCCCG", "TCCCT", "TCCGA", "TCCGC", "TCCGG", "TCCGT", "TCCTA", "TCCTC", "TCCTG",
	    "TCCTT", "TCGAA", "TCGAC", "TCGAG", "TCGAT", "TCGCA", "TCGCC", "TCGCG", "TCGCT", "TCGGA",
	    "TCGGC", "TCGGG", "TCGGT", "TCGTA", "TCGTC", "TCGTG", "TCGTT", "TCTAA", "TCTAC", "TCTAG",
	    "TCTAT", "TCTCA", "TCTCC", "TCTCG", "TCTCT", "TCTGA", "TCTGC", "TCTGG", "TCTGT", "TCTTA",
	    "TCTTC", "TCTTG", "TCTTT", "TGAAA", "TGAAC", "TGAAG", "TGAAT", "TGACA", "TGACC", "TGACG",
	    "TGACT", "TGAGA", "TGAGC", "TGAGG", "TGAGT", "TGATA", "TGATC", "TGATG", "TGATT", "TGCAA",
	    "TGCAC", "TGCAG", "TGCAT", "TGCCA", "TGCCC", "TGCCG", "TGCCT", "TGCGA", "TGCGC", "TGCGG",
	    "TGCGT", "TGCTA", "TGCTC", "TGCTG", "TGCTT", "TGGAA", "TGGAC", "TGGAG", "TGGAT", "TGGCA",
	    "TGGCC", "TGGCG", "TGGCT", "TGGGA", "TGGGC", "TGGGG", "TGGGT", "TGGTA", "TGGTC", "TGGTG",
	    "TGGTT", "TGTAA", "TGTAC", "TGTAG", "TGTAT", "TGTCA", "TGTCC", "TGTCG", "TGTCT", "TGTGA",
	    "TGTGC", "TGTGG", "TGTGT", "TGTTA", "TGTTC", "TGTTG", "TGTTT", "TTAAA", "TTAAC", "TTAAG",
	    "TTAAT", "TTACA", "TTACC", "TTACG", "TTACT", "TTAGA", "TTAGC", "TTAGG", "TTAGT", "TTATA",
	    "TTATC", "TTATG", "TTATT", "TTCAA", "TTCAC", "TTCAG", "TTCAT", "TTCCA", "TTCCC", "TTCCG",
	    "TTCCT", "TTCGA", "TTCGC", "TTCGG", "TTCGT", "TTCTA", "TTCTC", "TTCTG", "TTCTT", "TTGAA",
	    "TTGAC", "TTGAG", "TTGAT", "TTGCA", "TTGCC", "TTGCG", "TTGCT", "TTGGA", "TTGGC", "TTGGG",
	    "TTGGT", "TTGTA", "TTGTC", "TTGTG", "TTGTT", "TTTAA", "TTTAC", "TTTAG", "TTTAT", "TTTCA",
	    "TTTCC", "TTTCG", "TTTCT", "TTTGA", "TTTGC", "TTTGG", "TTTGT", "TTTTA", "TTTTC", "TTTTG",
	    "TTTTT" } }
};

/* Checks that the filepath is readable and exits if it is not. */
static inline void
assert_readable(const std::string& path)
{
	if (access(path.c_str(), R_OK) == -1) {
		std::cerr << PROGRAM ": error: `" << path << "': " << strerror(errno) << std::endl;
		exit(EXIT_FAILURE);
	}
}

/* Checks if the base is ATGC. */
bool
isAcceptedBase(unsigned char C)
{
	return (C == 'A' || C == 'T' || C == 'G' || C == 'C');
}

char
RC(unsigned char C)
{
	switch (C) {
	case 'A':
	case 'a':
		return 'T';
	case 'T':
	case 't':
		return 'A';
	case 'G':
	case 'g':
		return 'C';
	case 'C':
	case 'c':
		return 'G';
	default:
		return 'N';
	}
}

/* Find the first only ATGC kmer starting at the beginning of a sequence.
 * 	Assumption: no insertions or deletions. */
unsigned
findFirstAcceptedKmer(unsigned b_i, const std::string& contigSeq)
{
	for (unsigned i = b_i; i + opt::k < contigSeq.size();) {
		if (isAcceptedBase(toupper(contigSeq.at(i)))) {
			bool good_kmer = true;
			for (unsigned j = i + 1; j < i + opt::k; j++) {
				if (!isAcceptedBase(toupper(contigSeq.at(j)))) {
					good_kmer = false;
					i = j + 1;
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
	return contigSeq.size() - 1;
}

/* Helper for filling out the LPS array for detecting a low complexity repeat. */
void
computeLPSArray(std::string possible_repeat, int n, std::vector<int>& lps)
{
	int len = 0;
	int i;

	lps[0] = 0;
	i = 1;
	while (i < n) {
		if (possible_repeat[i] == possible_repeat[len]) {
			len++;
			lps[i] = len;
			i++;
		} else {
			if (len != 0) {
				len = lps[len - 1];
			} else {
				lps[i] = 0;
				i++;
			}
		}
	}
}

/* Determines if a string is a low complexity repeat of a word. */
bool
isRepeatInsertion(const std::string& possible_repeat)
{
	int n = possible_repeat.size();
	std::vector<int> lps(n);

	computeLPSArray(possible_repeat, n, lps);

	int len = lps[n - 1];
	return (len > 0 && n % (n - len) == 0);
}

/* Struct that keeps track of details for substitutions. */
struct sRec
{
	unsigned pos = 0;
	unsigned char draft_char = 0;
	unsigned char sub_base = 0;
	unsigned num_support = 0;
};

struct seqNode
{
	int node_type = -1; // -1=unset; 0=position; 1=character
	size_t s_pos = 0;
	size_t e_pos = 0;
	unsigned char c = 0;
	unsigned num_support = 0;
};

/* Makes a character insertion RIGHT BEFORE <insert_pos> by creating a character node holding <c>
 * with <num_support>
 * 	- sets <node> to your insertion node (the node that holds the character <c>). */
void
makeInsertion(
    unsigned& t_node_index,
    int insert_pos,
    const std::string& insertion_bases,
    unsigned num_support,
    std::vector<seqNode>& newSeq)
{
	seqNode orig_node = newSeq[t_node_index];
	std::vector<seqNode> to_insert;
	for (char insertion_base : insertion_bases) {
		seqNode insertion_node;
		insertion_node.node_type = 1;
		insertion_node.c = insertion_base;
		insertion_node.num_support = num_support;
		to_insert.push_back(insertion_node);
	}
	if (orig_node.node_type == 0) {
		if (insert_pos <= orig_node.s_pos) {
			// gather nodes following this insertion
			std::vector<seqNode> reappend;
			unsigned i = t_node_index;
			while (i < newSeq.size() && newSeq[i].node_type != -1) {
				reappend.push_back(newSeq[i]);
				newSeq[i].node_type = -1;
				i++;
			}
			// make insertion
			for (unsigned i = 0; i < to_insert.size(); i++) {
				if (t_node_index + i < newSeq.size()) {
					newSeq[t_node_index + i] = to_insert[i];
				} else {
					newSeq.push_back(to_insert[i]);
				}
			}
			// reappend
			for (unsigned i = 0; i < reappend.size(); i++) {
				if (t_node_index + to_insert.size() + i < newSeq.size()) {
					newSeq[t_node_index + to_insert.size() + i] = reappend[i];
				} else {
					newSeq.push_back(reappend[i]);
				}
			}
		} else {
			seqNode after_node;
			after_node.node_type = 0;
			after_node.s_pos = insert_pos;
			after_node.e_pos = orig_node.e_pos;
			newSeq[t_node_index].e_pos = insert_pos - 1;
			for (unsigned i = 0; i < to_insert.size(); i++) {
				if (t_node_index + i + 1 < newSeq.size()) {
					newSeq[t_node_index + i + 1] = to_insert[i];
				} else {
					newSeq.push_back(to_insert[i]);
				}
			}
			if (t_node_index + to_insert.size() + 1 < newSeq.size()) {
				newSeq[t_node_index + to_insert.size() + 1] = after_node;
			} else {
				newSeq.push_back(after_node);
			}
			t_node_index++;
		}
	} else if (orig_node.node_type == 1) {
		// gather the nodes following this insertion
		unsigned i = t_node_index;
		std::vector<seqNode> reappend;
		while (i < newSeq.size() && newSeq[i].node_type != -1) {
			reappend.push_back(newSeq[i]);
			newSeq[i].node_type = -1;
			i++;
		}
		// make the insertion
		for (unsigned i = 0; i < to_insert.size(); i++) {
			if (t_node_index + i < newSeq.size()) {
				newSeq[t_node_index + i] = to_insert[i];
			} else {
				newSeq.push_back(to_insert[i]);
			}
		}
		// push all of the insertions back
		for (unsigned i = 0; i < reappend.size(); i++) {
			if (t_node_index + to_insert.size() + i < newSeq.size()) {
				newSeq[t_node_index + to_insert.size() + i] = reappend[i];
			} else {
				newSeq.push_back(reappend[i]);
			}
		}
	}
}

/* Make a deletion starting and including <pos> of length <num_del> with support <num_support> in
 * seqNode structure.
 * 	- sets the t_node_index and pos to the position right after the deletion */
void
makeDeletion(
    unsigned& t_node_index,
    unsigned& pos,
    unsigned num_del,
    unsigned num_support,
    std::vector<seqNode>& newSeq)
{
	seqNode orig_node = newSeq[t_node_index];
	if (orig_node.node_type == 0) {
		unsigned leftover_del = 0;
		if (pos <= orig_node.s_pos) {
			if (pos + num_del <= orig_node.e_pos) {
				// we are deleting off the beginning of a position node
				newSeq[t_node_index].s_pos = pos + num_del;
				newSeq[t_node_index].num_support = num_support;
				pos = newSeq[t_node_index].s_pos;
				return;
			}
			// we deleted the entire position node and are moving on
			leftover_del = pos + num_del - orig_node.e_pos;
			pos = orig_node.e_pos + 1;
			// overwite the following onto this one
			unsigned i = t_node_index + 1;
			while (i < newSeq.size() && newSeq[i].node_type != -1) {
				newSeq[i - 1] = newSeq[i];
				newSeq[i].node_type = -1;
				i++;
			}
		} else {
			if (pos + num_del <= orig_node.e_pos) {
				// we are deleting in the middle of a position node
				seqNode split_node;
				split_node.node_type = 0;
				split_node.s_pos = pos + num_del;
				split_node.e_pos = orig_node.e_pos;
				split_node.num_support = num_support;
				newSeq[t_node_index].e_pos = pos - 1;
				pos = split_node.s_pos;
				t_node_index++;
				if (t_node_index < newSeq.size()) {
					newSeq[t_node_index] = split_node;
				} else {
					newSeq.push_back(split_node);
				}
				return;
			}
			// deleted from the middle of a position node past the end of it
			leftover_del = pos + num_del - orig_node.e_pos;
			newSeq[t_node_index].e_pos = pos - 1;
			pos = orig_node.e_pos + 1;
			t_node_index++;
		}
		if (leftover_del > 0) {
			// pass the deletion to the next seqNode
			if (t_node_index < newSeq.size() && newSeq[t_node_index].node_type != -1) {
				if (newSeq[t_node_index].node_type == 0) {
					pos = newSeq[t_node_index].s_pos;
				}
				makeDeletion(t_node_index, pos, leftover_del, num_support, newSeq);
			}
		}
	} else if (orig_node.node_type == 1) {
		unsigned i = t_node_index;
		unsigned leftover_del = num_del;
		// delete all the characters as possible
		while (i < newSeq.size() && newSeq[i].node_type == 1 && leftover_del > 0) {
			newSeq[i].node_type = -1;
			leftover_del--;
			i++;
		}
		// overwrite what comes after the characters
		unsigned j = t_node_index;
		while (i < newSeq.size() && newSeq[i].node_type != -1) {
			newSeq[j] = newSeq[i];
			newSeq[i].node_type = -1;
			i++;
			j++;
		}
		// deal with whatever is left
		if (leftover_del > 0) {
			// pass the deletion to the next seqNode
			if (t_node_index < newSeq.size() && newSeq[t_node_index].node_type != -1) {
				if (newSeq[t_node_index].node_type == 0) {
					pos = newSeq[t_node_index].s_pos;
				}
				makeDeletion(t_node_index, pos, leftover_del, num_support, newSeq);
			}
		}
	}
}

/* Returns the character at pos based on the seqNode structure. */
unsigned char
getCharacter(unsigned& pos, seqNode node, const string& contigSeq)
{
	if (node.node_type == 0) {
		return contigSeq.at(pos);
	}
	if (node.node_type == 1) {
		return node.c;
	}
	unsigned char c = 0;
	return c;
}

/* Increments the position and adjusts the node accordingly based on seqNode structure. */
void
increment(unsigned& pos, unsigned& node_index, vector<seqNode>& newSeq)
{
	seqNode node = newSeq[node_index];
	if (node.node_type == 0) {
		pos++;
		if (pos > node.e_pos) {
			node_index++;
			if (node_index < newSeq.size() && newSeq[node_index].node_type == 0) {
				pos = newSeq[node_index].s_pos;
			}
		}
	} else if (node.node_type == 1) {
		node_index++;
		if (node_index < newSeq.size() && newSeq[node_index].node_type == 0) {
			pos = newSeq[node_index].s_pos;
		}
	}
}

/* Find the first accepted kmer (contains only ATGC characters) starting anywhere, based on the
 * RopeLink structure. */
std::string
findAcceptedKmer(
    unsigned& h_seq_i,
    unsigned& t_seq_i,
    unsigned& h_node_index,
    unsigned& t_node_index,
    const string& contigSeq,
    vector<seqNode>& newSeq)
{
	// temporary values
	std::string kmer_str;
	seqNode curr_node = newSeq[t_node_index];
	unsigned temp_t_node_index = t_node_index;
	unsigned temp_h_node_index;
	unsigned i = t_seq_i;
	while (i < contigSeq.size() && temp_t_node_index < newSeq.size() &&
	       newSeq[temp_t_node_index].node_type != -1) {
		unsigned char c;
		c = getCharacter(i, curr_node, contigSeq);
		if (isAcceptedBase(toupper(c))) {
			std::string kmer_str;
			kmer_str += c;
			temp_h_node_index = temp_t_node_index;
			unsigned j = i;
			increment(j, temp_t_node_index, newSeq);
			// continue until you cant
			while (j < contigSeq.size() && temp_t_node_index < newSeq.size() &&
			       newSeq[temp_t_node_index].node_type != -1) {
				curr_node = newSeq[temp_t_node_index];
				c = getCharacter(j, curr_node, contigSeq);
				if (!isAcceptedBase(toupper(c))) {
					i = j;
					break;
				}
				kmer_str += c;
				if (kmer_str.size() == opt::k) {
					break;
				}
				increment(j, temp_t_node_index, newSeq);
			}
			// you found a good kmer so return it and adjust
			if (kmer_str.size() == opt::k) {
				h_seq_i = i;
				t_seq_i = j;
				h_node_index = temp_h_node_index;
				t_node_index = temp_t_node_index;
				return kmer_str;
			}
		}
		increment(i, temp_t_node_index, newSeq);
	}

	h_seq_i = contigSeq.length();
	t_seq_i = contigSeq.length();
	return "";
}

/* Get the previous insertion (aka continuous string of character nodes) starting at t_node_index.
 */
std::string
getPrevInsertion(unsigned t_seq_i, unsigned t_node_index, vector<seqNode>& newSeq)
{
	std::string prev_insertion;
	// if we just finished the insertion
	if ((t_node_index < newSeq.size() && newSeq[t_node_index].node_type == 0 &&
	     t_seq_i == newSeq[t_node_index].s_pos) ||
	    newSeq[t_node_index].node_type == 1) {
		t_node_index--;
	}
	while (t_node_index < newSeq.size() && newSeq[t_node_index].node_type == 1) {
		prev_insertion += RC(newSeq[t_node_index].c);
		t_node_index--;
	}
	return prev_insertion;
}

/* Write the edits and new draft contig into respective files. */
void
writeEditsToFile(
    std::ofstream& dfout,
    std::ofstream& rfout,
    const std::string& contigHdr,
    const std::string& contigSeq,
    std::vector<seqNode>& newSeq,
    std::queue<sRec>& substitution_record)
{
	dfout << ">" << contigHdr.c_str() << "\n";
	unsigned node_index = 0;
	std::string insertion_bases;
	int num_support = -1;
	unsigned char draft_char;
	unsigned pos = 0;

	// track a deletion
	std::string deleted_bases;
	seqNode curr_node = newSeq[node_index];
	while (node_index < newSeq.size() && curr_node.node_type != -1) {
		if (curr_node.node_type == 0) {
			draft_char = contigSeq.at(curr_node.s_pos);
			// log an insertion if it occured before this
			if (!insertion_bases.empty()) {
				rfout << contigHdr.c_str() << "\t" << pos + 1 << "\t" << draft_char << "\t+"
				      << insertion_bases.c_str() << "\t" << num_support << "\n";
				insertion_bases = "";
				num_support = -1;
			}
			// log all the substitutions up to this point
			while (!substitution_record.empty() &&
			       substitution_record.front().pos <= curr_node.e_pos) {
				rfout << contigHdr.c_str() << "\t" << substitution_record.front().pos + 1 << "\t"
				      << substitution_record.front().draft_char << "\t"
				      << substitution_record.front().sub_base << "\t"
				      << substitution_record.front().num_support << "\n";
				substitution_record.pop();
			}
			// fprintf(dfout, "%s", contigSeq.substr(curr_node.s_pos,
			// (curr_node.e_pos-curr_node.s_pos+1)).c_str());
			dfout << contigSeq.substr(curr_node.s_pos, (curr_node.e_pos - curr_node.s_pos + 1))
			             .c_str();
			pos = curr_node.e_pos + 1;
		} else if (curr_node.node_type == 1) {
			insertion_bases += curr_node.c;
			if (num_support == -1) {
				num_support = curr_node.num_support;
			}
			dfout << curr_node.c;
		}
		node_index++;
		if (node_index < newSeq.size()) {
			curr_node = newSeq[node_index];
			if (curr_node.node_type == 0 && curr_node.s_pos != pos) {
				// print out the deletion
				rfout << contigHdr.c_str() << "\t" << pos + 1 << "\t" << contigSeq.at(pos) << "\t-"
				      << contigSeq.substr(pos, (curr_node.s_pos - pos)).c_str() << "\t"
				      << curr_node.num_support << "\n";
			}
		}
	}
	dfout << "\n";
}

/* Roll ntHash using the seqNode structure. */
bool
roll(
    unsigned& h_seq_i,
    unsigned& t_seq_i,
    unsigned& h_node_index,
    unsigned& t_node_index,
    const std::string& contigSeq,
    std::vector<seqNode>& newSeq,
    unsigned char& charOut,
    unsigned char& charIn)
{

	// quit if h_seq_i is out of scope
	if (h_seq_i >= contigSeq.size() || h_node_index >= newSeq.size()) {
		return false;
	}
	charOut = getCharacter(h_seq_i, newSeq[h_node_index], contigSeq);
	increment(h_seq_i, h_node_index, newSeq);

	increment(t_seq_i, t_node_index, newSeq);
	// quit if t_seq_i is out of scope
	if (t_seq_i >= contigSeq.size() || t_node_index >= newSeq.size()) {
		return false;
	}
	charIn = getCharacter(t_seq_i, newSeq[t_node_index], contigSeq);

	return true;
}

/* Accept the edit */
void
makeEdit(
    unsigned char& draft_char,
    unsigned& best_edit_type,
    unsigned char& best_sub_base,
    string& best_indel,
    unsigned& best_num_support,
    std::queue<sRec>& substitution_record,
    unsigned& h_seq_i,
    unsigned& t_seq_i,
    unsigned& h_node_index,
    unsigned& t_node_index,
    uint64_t& fhVal,
    uint64_t& rhVal,
    uint64_t* hVal,
    std::string& contigSeq,
    std::vector<seqNode>& newSeq)
{

	bool skipped_repeat = false;
	std::string prev_insertion;
	// make our edit
	seqNode tNode = newSeq[t_node_index];
	switch (best_edit_type) {
	case 1: // SUBSTITUTION MADE
		// apply the change to the actual contigSequence
		if (tNode.node_type == 0) {
			contigSeq[t_seq_i] = best_sub_base;
			sRec subst;
			subst.draft_char = draft_char;
			subst.pos = t_seq_i;
			subst.sub_base = best_sub_base;
			subst.num_support = best_num_support;
			substitution_record.push(subst);
		} else if (tNode.node_type == 1) {
			newSeq[t_node_index].c = best_sub_base;
		}
		// make sure we change our current hash to match it
		NTMC64_changelast(draft_char, best_sub_base, opt::k, opt::h, fhVal, rhVal, hVal);
		if (opt::verbose) {
			std::cout << "\tt_seq_i: " << t_seq_i << " SUB: " << best_sub_base
			          << " check_present: " << best_num_support << std::endl;
		}
		break;
	case 2: // INSERTION MADE
		// check if we need to check the insertion
		// 	or low complexity and make that check before preceding
		prev_insertion = getPrevInsertion(t_seq_i, t_node_index, newSeq);
		if (prev_insertion.size() + best_indel.size() >= opt::k) {
			// check if the original previous insertion was a low complexity repeat or we have
			// reached our hard cap
			if (isRepeatInsertion(prev_insertion) ||
			    prev_insertion.size() + best_indel.size() >= opt::insertion_cap) {
				unsigned j = 1;
				if (newSeq[t_node_index].node_type == 0 && t_seq_i == newSeq[t_node_index].s_pos) {
					j = 0;
				}
				for (unsigned i = prev_insertion.size(); i > 0; i--) {
					if (t_node_index + j < newSeq.size() &&
					    newSeq[t_node_index + j].node_type != -1) {
						newSeq[t_node_index - i] = newSeq[t_node_index + j];
						newSeq[t_node_index + j].node_type = -1;
						j++;
					} else {
						newSeq[t_node_index - i].node_type = -1;
					}
				}
				NTMC64(
				    findAcceptedKmer(
				        h_seq_i, t_seq_i, h_node_index, t_node_index, contigSeq, newSeq)
				        .c_str(),
				    opt::k,
				    opt::h,
				    fhVal,
				    rhVal,
				    hVal);
				skipped_repeat = true;
			} else {
				// check the rest of the insertion with the extra insertion base(s) for low
				// complexity
				for (unsigned w = 0; w < best_indel.size();
				     w++) { // NOLINT(bugprone-too-small-loop-variable)
					prev_insertion.insert(prev_insertion.begin(), RC(best_indel[w]));
					if (isRepeatInsertion(prev_insertion)) {
						unsigned j = 1;
						if (newSeq[t_node_index].node_type == 0 &&
						    t_seq_i == newSeq[t_node_index].s_pos) {
							j = 0;
						}
						for (unsigned i = prev_insertion.size() - w; i > 0; i--) {
							if (t_node_index + j < newSeq.size() &&
							    newSeq[t_node_index + j].node_type != -1) {
								newSeq[t_node_index - i] = newSeq[t_node_index + j];
								newSeq[t_node_index + j].node_type = -1;
								j++;
							} else {
								newSeq[t_node_index - i].node_type = -1;
							}
						}
						NTMC64(
						    findAcceptedKmer(
						        h_seq_i, t_seq_i, h_node_index, t_node_index, contigSeq, newSeq)
						        .c_str(),
						    opt::k,
						    opt::h,
						    fhVal,
						    rhVal,
						    hVal);
						skipped_repeat = true;
					}
				}
			}
		}
		// if we didn't skip this region for a repeat, then make this insertion
		if (!skipped_repeat) {
			makeInsertion(t_node_index, t_seq_i, best_indel, best_num_support, newSeq);
			NTMC64_changelast(draft_char, best_indel[0], opt::k, opt::h, fhVal, rhVal, hVal);
			if (opt::verbose) {
				std::cout << "\tt_seq_i: "
				          << t_seq_i
				          //						<< " " << t_node_index
				          << " INS: " << best_indel << " check_present: " << best_num_support
				          << std::endl;
			}
		}
		break;
	case 3: // DELETION MADE
		if (opt::verbose) {
			std::cout << "\tt_seq_i: " << t_seq_i << " DEL: " << best_indel
			          << " check_present: " << best_num_support << std::endl;
		}
		makeDeletion(t_node_index, t_seq_i, best_indel.size(), best_num_support, newSeq);
		NTMC64_changelast(
		    draft_char,
		    getCharacter(t_seq_i, newSeq[t_node_index], contigSeq),
		    opt::k,
		    opt::h,
		    fhVal,
		    rhVal,
		    hVal);
		break;
	case 0:
		if (opt::verbose) {
			std::cout << "\tt_seq_i: " << t_seq_i << " FIX NOT FOUND" << std::endl;
		}
		break;
	default:
		break;
	}
}

/* Try a deletion in ntEdit. */
int
tryDeletion(
    const unsigned char draft_char,
    unsigned num_deletions,
    unsigned& h_seq_i,
    unsigned& t_seq_i,
    unsigned& h_node_index,
    unsigned& t_node_index,
    uint64_t& fhVal,
    uint64_t& rhVal,
    uint64_t* hVal,
    const std::string& contigSeq,
    std::vector<seqNode>& newSeq,
    BloomFilter& bloom,
    std::string& deleted_bases)
{

	// set temporary values
	uint64_t temp_fhVal = fhVal;
	uint64_t temp_rhVal = rhVal;
	unsigned temp_h_seq_i = h_seq_i;
	unsigned temp_t_seq_i = t_seq_i;
	unsigned temp_h_node_index = h_node_index;
	unsigned temp_t_node_index = t_node_index;
	unsigned char charOut;
	unsigned char charIn;

	// make the deletion
	for (unsigned i = 0; i < num_deletions; i++) {
		deleted_bases += getCharacter(temp_t_seq_i, newSeq[temp_t_node_index], contigSeq);
		increment(temp_t_seq_i, temp_t_node_index, newSeq);
	}
	NTMC64_changelast(
	    draft_char,
	    getCharacter(temp_t_seq_i, newSeq[temp_t_node_index], contigSeq),
	    opt::k,
	    opt::h,
	    temp_fhVal,
	    temp_rhVal,
	    hVal);

	// verify the deletion with a subset
	unsigned check_present = 0;
	if (bloom.contains(hVal)) {
		check_present++; // check for changing the kmer after deletion
	}
	for (unsigned k = 1; k <= (opt::k - 2) && temp_h_seq_i < contigSeq.size(); k++) {
		if (roll(
		        temp_h_seq_i,
		        temp_t_seq_i,
		        temp_h_node_index,
		        temp_t_node_index,
		        contigSeq,
		        newSeq,
		        charOut,
		        charIn)) {
			NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
			if (k%opt::jump == 0 && bloom.contains(hVal)) {
				check_present++;
			}
		}
	}

	if (opt::verbose) {
		std::cout << "\t\tdeleting: " << deleted_bases << " check_present: " << check_present
		          << std::endl;
	}
        if ((!opt::use_ratio && static_cast<float>(check_present) >= ( static_cast<float>(opt::k) / opt::edit_threshold))
		|| (opt::use_ratio && static_cast<float>(check_present) >= ( 1 + ( static_cast<float>(opt::k) / opt::jump )) * opt::edit_ratio)) { // RLW
		return static_cast<int>(check_present);
	}
	return 0;
}

/* Try indel combinations starting with index_char. */
bool
tryIndels(
    const unsigned char draft_char,
    const unsigned char index_char,
    unsigned& num_deletions,
    unsigned& h_seq_i,
    unsigned& t_seq_i,
    unsigned& h_node_index,
    unsigned& t_node_index,
    uint64_t& fhVal,
    uint64_t& rhVal,
    uint64_t* hVal,
    const std::string& contigSeq,
    std::vector<seqNode>& newSeq,
    BloomFilter& bloom,
    unsigned& best_edit_type,
    std::string& best_indel,
    unsigned& best_num_support)
{

	// initialize temporary values
	uint64_t temp_fhVal;
	uint64_t temp_rhVal;
	unsigned temp_h_seq_i;
	unsigned temp_t_seq_i;
	unsigned temp_h_node_index;
	unsigned temp_t_node_index;
	unsigned temp_best_num_support = 0;
	std::string temp_best_indel;
	unsigned temp_best_edit_type = 0;
	unsigned char charIn;
	unsigned char charOut;

	// try all of the combinations of indels starting with our index_char
	for (unsigned i = 0; i < num_tries[opt::max_insertions]; i++) {
		// gather the insertion bases
		std::string insertion_bases = multi_possible_bases[index_char][i];
		insertion_bases += draft_char;

		// set temporary values
		temp_fhVal = fhVal;
		temp_rhVal = rhVal;
		temp_h_seq_i = h_seq_i;
		temp_t_seq_i = t_seq_i;
		temp_h_node_index = h_node_index;
		temp_t_node_index = t_node_index;

		// change the last base
		NTMC64_changelast(draft_char, index_char, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
		unsigned check_present = 0;
		unsigned k = 0; // RLW
		// check subset with the insertion
		for (; k < insertion_bases.size() - 1 && temp_h_seq_i < contigSeq.size(); k++) {
			NTMC64(
			    getCharacter(temp_h_seq_i, newSeq[temp_h_node_index], contigSeq),
			    insertion_bases[k+1],
			    opt::k,
			    opt::h,
			    temp_fhVal,
			    temp_rhVal,
			    hVal);
			increment(temp_h_seq_i, temp_h_node_index, newSeq);
			if (k % opt::jump == 0 && bloom.contains(hVal)) {// RLW
				check_present++;
			}
		}
		// check subset after insertion
		for (; k < opt::k - 1 && temp_h_seq_i < contigSeq.size(); k++) {
			if (roll(
			        temp_h_seq_i,
			        temp_t_seq_i,
			        temp_h_node_index,
			        temp_t_node_index,
			        contigSeq,
			        newSeq,
			        charOut,
			        charIn)) {
				NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
				if (k % opt::jump == 0 && bloom.contains(hVal)) {//RLW
					check_present++;
				}
			}
		}
		insertion_bases.pop_back();
		if (opt::verbose) {
			std::cout << "\t\tinserting: " << insertion_bases << " check_present: " << check_present
			          << std::endl;
		}
		// if the insertion is good, store the insertion accordingly RLW
		if ((!opt::use_ratio && static_cast<float>(check_present) >= (static_cast<float>(opt::k) / opt::edit_threshold)) || (opt::use_ratio && static_cast<float>(check_present) >= ( static_cast<float>(opt::k) / opt::jump) * opt::edit_ratio)) { // RLW
			if (opt::mode == 0) {
				// if we are in default mode, we just accept this first good insertion and return
				best_edit_type = 2;
				best_indel = insertion_bases;
				best_num_support = check_present;
				return true;
			}
			if (opt::mode == 1 || opt::mode == 2) {
				// if we are in some deep mode, we look for the best indel within index char first
				if (check_present > temp_best_num_support) {
					temp_best_edit_type = 2;
					temp_best_indel = insertion_bases;
					temp_best_num_support = check_present;
				}
			}
		}

		if (num_deletions <= opt::max_deletions) {
			std::string deleted_bases;
			unsigned del_support = tryDeletion(
			    draft_char,
			    num_deletions,
			    h_seq_i,
			    t_seq_i,
			    h_node_index,
			    t_node_index,
			    fhVal,
			    rhVal,
			    hVal,
			    contigSeq,
			    newSeq,
			    bloom,
			    deleted_bases);
			if (del_support > 0) {
				if (opt::mode == 0) {
					best_edit_type = 3;
					best_indel = deleted_bases;
					best_num_support = del_support;
					return true;
				}
				if (opt::mode == 1 || opt::mode == 2) {
					if (del_support > temp_best_num_support) {
						temp_best_edit_type = 3;
						temp_best_indel = deleted_bases;
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
void
kmerizeAndCorrect(
    string& contigHdr,
    string& contigSeq,
    unsigned seqLen,
    BloomFilter& bloom,
    std::ofstream& dfout,
    std::ofstream& rfout)
{

	// initialize values for hashing
	uint64_t fhVal = 0;
	uint64_t rhVal = 0;
	uint64_t* hVal;
	unsigned char charIn = 0;
	unsigned char charOut;
	unsigned char draft_char;
	hVal = new uint64_t[opt::h];

	// vector to record substitutions
	std::queue<sRec> substitution_record;

	// initialize and readjust the first character depending on the first N or nonATGC kmer
	unsigned h_seq_i = findFirstAcceptedKmer(0, contigSeq);
	unsigned t_seq_i = h_seq_i + opt::k - 1;

	// intialize our seed kmer
	if (h_seq_i + opt::k - 1 < seqLen) {
		NTMC64(contigSeq.substr(h_seq_i, opt::k).c_str(), opt::k, opt::h, fhVal, rhVal, hVal);
		charIn = contigSeq.at(t_seq_i);
	}

	// compressed string
	//	vector<seqNode> newSeq(contigSeq.size()/4);
	std::vector<seqNode> newSeq;
	newSeq.reserve(contigSeq.size() / 4);
	// initialize our first node
	seqNode root;
	root.node_type = 0;
	root.s_pos = 0;
	root.e_pos = seqLen - 1;
	//	newSeq[0] = root;
	newSeq.push_back(root);
	// h_seq_i and t_seq_i node pointers
	unsigned h_node_index = 0;
	unsigned t_node_index = 0;

	bool continue_edit = true;
	do {
		if (h_seq_i + opt::k - 1 >= seqLen) {
			break;
		}
		if (opt::verbose) {
			std::cout << h_seq_i << " " << t_seq_i << " " << charIn << " " << h_node_index << " "
			          << t_node_index << " " << hVal[0] << hVal[1] << hVal[2] << std::endl;
		}
		if (!bloom.contains(hVal)) {
			// make temporary value holders
			uint64_t temp_fhVal = fhVal;
			uint64_t temp_rhVal = rhVal;
			unsigned temp_h_seq_i = h_seq_i;
			unsigned temp_t_seq_i = t_seq_i;
			unsigned temp_h_node_index = h_node_index;
			unsigned temp_t_node_index = t_node_index;

			// set draft char if we choose to make an edit later
			draft_char = toupper(charIn);

			// confirm missing by checking subset
			unsigned check_missing = 0;
			bool do_not_fix = false;
			for (unsigned k = 0; k < opt::k && temp_h_seq_i < seqLen; k++) { // RLW
				if (roll(
				        temp_h_seq_i,
				        temp_t_seq_i,
				        temp_h_node_index,
				        temp_t_node_index,
				        contigSeq,
				        newSeq,
				        charOut,
				        charIn)) {
					NTMC64(charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
					if (!isAcceptedBase(toupper(charIn))) {
						do_not_fix = true;
						break;
					}
					if (k % opt::jump == 0 && !bloom.contains(hVal)) { // RLW
						check_missing++;
					}
				} else {
					do_not_fix = true;
					break;
				}
			}

			if (opt::verbose) {
				std::cout << "\tcheck_missing: " << check_missing << std::endl;
			}
			if (!do_not_fix && ((!opt::use_ratio && static_cast<float>(check_missing) >= (static_cast<float>(opt::k) / opt::missing_threshold))
				|| (opt::use_ratio && static_cast<float>(check_missing) >= (( static_cast<float>(opt::k) / opt::jump) * opt::missing_ratio)))) { // RLW
				// recorders
				unsigned num_deletions = 1;
				unsigned best_edit_type =
				    0; // 0=no edit made; 1=substitution; 2=insertion; 3=deletions
				std::string best_indel;
				unsigned char best_sub_base;
				unsigned best_num_support = 0;

				// try substitution
				for (const unsigned char sub_base : bases_array[draft_char]) {
					// reset the temporary values
					temp_fhVal = fhVal;
					temp_rhVal = rhVal;

					// hash the substitution change
					NTMC64_changelast(
					    draft_char, sub_base, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);

					// only do verification of substition if it is found in bloom filter
					if (bloom.contains(hVal) || opt::mode == 2) {
						// reset temporary values
						temp_h_node_index = h_node_index;
						temp_t_node_index = t_node_index;
						temp_h_seq_i = h_seq_i;
						temp_t_seq_i = t_seq_i;

						// change the substitution
						if (newSeq[t_node_index].node_type == 0) {
							contigSeq.at(temp_t_seq_i) = sub_base;
						} else if (newSeq[t_node_index].node_type == 1) {
							newSeq[t_node_index].c = sub_base;
						}
						// check the subset to see if this substition is good
						unsigned check_present = 0;
						for (unsigned k = 0; // RLW
						     k < opt::k && temp_h_seq_i < seqLen && temp_t_seq_i < seqLen;
						     k++) {
							if (roll(
							        temp_h_seq_i,
							        temp_t_seq_i,
							        temp_h_node_index,
							        temp_t_node_index,
							        contigSeq,
							        newSeq,
							        charOut,
							        charIn)) {
								NTMC64(
								    charOut, charIn, opt::k, opt::h, temp_fhVal, temp_rhVal, hVal);
								if (k % opt::jump == 0 && bloom.contains(hVal)) { // RLW
									check_present++;
								}
							} else {
								break;
							}
						}

						// revert the substitution
						if (newSeq[t_node_index].node_type == 0) {
							contigSeq.at(t_seq_i) = draft_char;
						} else if (newSeq[t_node_index].node_type == 1) {
							newSeq[t_node_index].c = draft_char;
						}
						if (opt::verbose) {
							std::cout << "\t\tsub: " << sub_base
							          << " check_present: " << check_present << std::endl;
						}

						if ((!opt::use_ratio && static_cast<float>(check_present) >= (static_cast<float>(opt::k) / opt::edit_threshold))
							|| (opt::use_ratio && static_cast<float>(check_present) >= (( static_cast<float>(opt::k) / opt::jump))*opt::edit_ratio)) { // RLW

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
							if (tryIndels(
							        draft_char,
							        sub_base,
							        num_deletions,
							        h_seq_i,
							        t_seq_i,
							        h_node_index,
							        t_node_index,
							        fhVal,
							        rhVal,
							        hVal,
							        contigSeq,
							        newSeq,
							        bloom,
							        best_edit_type,
							        best_indel,
							        best_num_support)) {
								if (opt::mode == 0 || opt::mode == 1) {
									break;
								}
							}
						}
					}
				}

				makeEdit(
				    draft_char,
				    best_edit_type,
				    best_sub_base,
				    best_indel,
				    best_num_support,
				    substitution_record,
				    h_seq_i,
				    t_seq_i,
				    h_node_index,
				    t_node_index,
				    fhVal,
				    rhVal,
				    hVal,
				    contigSeq,
				    newSeq);
			}
		}
		// roll and skip over non-ATGC containing kmers
		int target_t_seq_i = -1;
		do {
			if (roll(
			        h_seq_i,
			        t_seq_i,
			        h_node_index,
			        t_node_index,
			        contigSeq,
			        newSeq,
			        charOut,
			        charIn)) {
				if (!isAcceptedBase(toupper(charIn))) {
					target_t_seq_i = static_cast<int>(t_seq_i) + static_cast<int>(opt::k);
				}
				NTMC64(charOut, charIn, opt::k, opt::h, fhVal, rhVal, hVal);
			} else {
				continue_edit = false;
				break;
			}
		} while (target_t_seq_i >= 0 && t_seq_i != target_t_seq_i);
	} while (continue_edit);

	// clean allocated memory for hash
	delete[] hVal;
	hVal = nullptr;

#pragma omp critical(write)
	{
		// write this to file
		writeEditsToFile(dfout, rfout, contigHdr, contigSeq, newSeq, substitution_record);
	}
}

/* Read the contigs from the file and polish each contig. */
void
readAndCorrect(BloomFilter& bloom)
{
	// read file handle
	gzFile dfp;
	dfp = gzopen(opt::draft_filename.c_str(), "r");
	kseq_t* seq = kseq_init(dfp);
	int num_contigs = 0;
	constexpr int print_step_size = 1000000;

	// outfile handles
	std::string d_filename = opt::outfile_prefix + "_edited.fa";
	std::string r_filename = opt::outfile_prefix + "_changes.tsv";
	ofstream dfout;
	ofstream rfout;
	dfout.open(d_filename);
	rfout.open(r_filename);

	rfout << "ID\tbpPosition+1\tOriginalBase\tNewBase\tSupport " << opt::k << "-mer (out of "
	      << ((opt::k / opt::jump) + 1) << ")\tAlternateNewBase\tAlt.Support " << opt::k << "-mers\n"; // RLW

#pragma omp parallel shared(seq, dfout, rfout)
	{
		std::string contigHdr;
		std::string contigSeq;
		std::string contigName;
		bool stop = false;

		while (true) {
#pragma omp critical(reading)
			{
				if (!stop && kseq_read(seq) >= 0) {
					contigHdr = seq->name.s;
					if (seq->comment.l) {
						contigName = contigHdr + " " + seq->comment.s;
					} else {
						contigName = contigHdr;
					}
					contigSeq = seq->seq.s;
				} else {
					stop = true;
				}
			}
			if (stop) {
				break;
			}
			unsigned seq_len = contigSeq.length();
			if (opt::verbose) {
				std::cout << contigName << std::endl;
			}
			if (seq_len >= opt::min_contig_len) {
				kmerizeAndCorrect(contigName, contigSeq, seq_len, bloom, dfout, rfout);
			}
#pragma omp atomic
			num_contigs++;
			if (num_contigs % print_step_size == 0) {
				std::cout << "Processed " << num_contigs << std::endl;
			}
		}
	}
	//#pragma omp barrier
	kseq_destroy(seq);
	gzclose(dfp);
	dfout.close();
	rfout.close();
}

int
main(int argc, char** argv)
{
	bool die = false;
	for (int c; (c = getopt_long(argc, argv, shortopts, longopts, nullptr)) != -1;) {
		std::istringstream arg(optarg != nullptr ? optarg : "");
		switch (c) {
		case '?':
			die = true;
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
		case 'X':
			arg >> opt::missing_ratio;
			opt::use_ratio=true;
			break;
		case 'Y':
			arg >> opt::edit_ratio;
			opt::use_ratio=true;
			break;
		case 'c':
			arg >> opt::insertion_cap;
			break;
		case 'j':
			arg >> opt::jump;
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
		default:
			break;
		}
		if (optarg != nullptr && (!arg.eof() || arg.fail())) {
			std::cerr << PROGRAM ": invalid option: `-" << static_cast<char>(c) << optarg << "'\n";
			exit(EXIT_FAILURE);
		}
	}

	// std::cout << opt::nthreads << " = thread no" << std::endl;

	time_t rawtime;
	time(&rawtime);
	std::cout << "\n-----------------Running ntEdit------------- " << ctime(&rawtime);

	// check the draft file is specified
	if (opt::draft_filename.empty()) {
		std::cerr << PROGRAM ": error: need to specify assembly draft file (-f)\n";
		die = true;
	} else {
		// if the file is specified check that it is readable
		assert_readable(opt::draft_filename);
	}

	// check that the bloom filter file is specified
	if (opt::bloom_filename.empty()) {
		std::cerr << PROGRAM ": error: need to specify the bloom filter file (-r)\n";
		die = true;
	} else {
		// if the file is specified check that it is readable
		assert_readable(opt::bloom_filename);
	}

	if (die) {
		std::cerr << "Try `" << PROGRAM << " --help' for more information.\n";
		exit(EXIT_FAILURE);
	}

	// check that the parameters x and y are in bound
	if (opt::missing_threshold < 3 && opt::missing_threshold > static_cast<float>(opt::k) &&
	    opt::edit_threshold < 3 && opt::edit_threshold > static_cast<float>(opt::k)) {
		std::cerr << PROGRAM ": warning: x and y parameters must be >=3 and <=k; x and y were "
		                     "reset to default values x=5, y=9.\n";
		constexpr float missing_threshold_override = 5;
		constexpr float edit_threshold_override = 5;
		opt::missing_threshold = missing_threshold_override;
		opt::edit_threshold = edit_threshold_override;
	}

	// check that the parameters i and d are in bound
	if ((opt::max_insertions == 0 && opt::max_deletions > 0) ||
	    (opt::max_insertions == 1 && opt::max_deletions > 1)) {
		std::cerr << PROGRAM ": warning: i and d parameter combination is not possible; d was set "
		                     "to the value of i.\n";
		opt::max_deletions = opt::max_insertions;
	}

	// get the basename for the file
	std::string draft_basename =
	    opt::draft_filename.substr(opt::draft_filename.find_last_of("/\\") + 1);
	std::string bloom_basename =
	    opt::bloom_filename.substr(opt::bloom_filename.find_last_of("/\\") + 1);

	// set the outfile prefix if it wasn't given
	if (opt::outfile_prefix.empty()) {
		std::ostringstream outfile_name;
		outfile_name << draft_basename << "_k" << opt::k << "_z" << opt::min_contig_len << "_r"
		             << bloom_basename << "_i" << opt::max_insertions << "_d" << opt::max_deletions
		             << "_m" << opt::mode;
		opt::outfile_prefix = outfile_name.str();
	}

	// print parameters:
        if (opt::use_ratio) { // RLW
                std::cout << "Running: " << PROGRAM << "\n -t " << opt::nthreads << "\n -f " << draft_basename
                << "\n -k " << opt::k << "\n -z " << opt::min_contig_len << "\n -b "
                << opt::outfile_prefix << "\n -r " << bloom_basename << "\n -i "
                << opt::max_insertions << "\n -d " << opt::max_deletions << "\n -X "
                << opt::missing_ratio << "\n -Y " << opt::edit_ratio << "\n -j " << opt::jump << "\n -m " << opt::mode
                << "\n -v " << opt::verbose << std::endl;
        } else {
		std::cout << "Running: " << PROGRAM << "\n -t " << opt::nthreads << "\n -f " << draft_basename
		<< "\n -k " << opt::k << "\n -z " << opt::min_contig_len << "\n -b "
		<< opt::outfile_prefix << "\n -r " << bloom_basename << "\n -i "
		<< opt::max_insertions << "\n -d " << opt::max_deletions << "\n -x "
		<< opt::missing_threshold << "\n -y " << opt::edit_threshold << "\n -j " << opt::jump << "\n -m " << opt::mode
		<< "\n -v " << opt::verbose << std::endl;
	}

	// Threading information
	omp_set_num_threads(static_cast<int>(opt::nthreads));

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
	std::cout << "\n----------Reading and Polishing Draft---------- " << ctime(&rawtime);
	readAndCorrect(bloom);

	time(&rawtime);
	std::cout << "\n----------ntEdit Polishing Complete !---------- " << ctime(&rawtime);

	return 0;
}
