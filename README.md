![Logo](https://github.com/bcgsc/ntEdit/blob/master/ntedit-logo.png)

# ntEdit

## Scalable genome sequence polishing
## 2/2019
## email: rwarren [at] bcgsc [dot] ca


### Description
-----------

ntEdit is a scalable genomics application for polishing genome assembly drafts.
ntEdit simplifies polishing and "haploidization" of gene and genome sequences with its re-usable Bloom filter design.
We expect ntEdit to have additional applications in fast mapping of simple nucleotide variations between any two individuals or species’ genomes.


### Implementation and requirements
-------------------------------

ntEdit v1.2.0 is written in C++. 

(We compiled with gcc 5.5.0)


### Install
-------

Clone this directory and enter this directory.
<pre>
https://github.com/bcgsc/ntEdit.git
cd ntEdit
</pre>
Compile ntEdit.
<pre>
make ntedit
</pre>


### Dependencies
-------

1. ntHits (https://github.com/bcgsc/nthits)
2. BloomFilter utilities (provided in ./lib)
3. kseq (provided in ./lib)


### Documentation
-------------

Refer to the README.md file on how to run ntEdit and our manuscript for information about the software and its performance 
Questions or comments?  We would love to hear from you!
rwarren@bcgsc.ca


### Citing ntEdit 
------------

<pre>
ntEdit: scalable genome sequence polishing
Rene L Warren, Lauren Coombe, Hamid Mohamadi, Jessica Zhang, Barry Jaquish, Nathalie Isabel, Steven JM Jones, Jean Bousquet, Joerg Bohlmann, Inanc Birol. 
2019. Bioinformatics (accepted issue TBD)
bioRxiv 2019; 565374; doi: https://doi.org/10.1101/565374

</pre>

The experimental data described in our paper can be downloaded here: http://www.bcgsc.ca/downloads/btl/ntedit/
Thank you for using, developing and promoting this free software.


### Credits
-------

ntEdit (concept, algorithm design and prototype):
Rene Warren

nthits / nthash / BloomFilter.pm:
Hamid Mohamadi

C++ implementation:
Jessica Zhang


### How to run in a pipeline
-------

1. Running nthits (please see nthits documentation : https://github.com/bcgsc/ntHits)
<pre>
nthits -c <kmer coverage threshold> -b <Bloom filter bit size> -k <kmer length> -t <number of threads> reads
eg.
nthits -c 2 --outbloom -p solidBF -b 36 -k 25 -t 48 Sim_100_300_1.fq Sim_100_300_2.fq
or
nthits -c 2 --outbloom -p solidBF -b 36 -k 25 -t 48 @reads.in

Where @reads.in is a file listing the path to all read fastq files to kmerize
note the options --outbloom and -p that must be set when using ntHits with ntEdit. 

If not specifying a hard threshold (-c), and relying instead on ntCard* to identify the error kmer coverage, please run with the --solid:

nthits -b 36 -k 50 -t 48 --outbloom --solid @reads.in 

NOTE: THIS WILL WORK WELL WITH ntEdit ONLY IF YOU HAVE SUFFICIENT READ COVERAGE (>30X), OTHERWISE SET KMER COVERAGE TO -c2 (>=20X) or -c 1 (<20X).

*Bioinformatics. 2017 May 1; 33(9): 1324–1330.
Published online 2017 Jan 5. doi: 10.1093/bioinformatics/btw832
PMCID: PMC5408799
PMID: 28453674
ntCard: a streaming algorithm for cardinality estimation in genomics data
Hamid Mohamadi, Hamza Khan and Inanc Birol
</pre>


2. Running ntEdit (see complete usage below)
<pre>
./ntedit -f <fasta file to polish> -k <kmer length> -r <Bloom filter from nthits> -b <base output name> -t <threads>

eg.
./ntedit -f ecoliWithMismatches001Indels0001.fa -r solidBF_k25.bf -k 25 -b ntEditEcolik25 -t 48
</pre>

### Running ntEdit
-------------
<pre>
e.g. ./ntedit -f ecoliWithMismatches001Indels0001.fa -r solidBF_k25.bf -k 25 -b ntEditEcolik25

 Options:
	-t,	number of threads [default=1]
	-f,	Draft genome assembly (FASTA, Multi-FASTA, and/or gzipped compatible), REQUIRED
	-r,	Bloom filter file (generated from ntHits), REQUIRED
	-b,	output file prefix, OPTIONAL
	-k,	kmer size, REQUIRED
	-z,	minimum contig length [default=100]
	-i,	maximum number of insertion bases to try, range 0-5, [default=4]
	-d,	maximum number of deletions bases to try, range 0-5, [default=5]
	-x,	k/x ratio for the number of kmers that should be missing, [default=5.000]
	-y, 	k/y ratio for the number of editted kmers that should be present, [default=9.000]
	-c,	cap for the number of base insertions that can be made at one position, [default=k*1.5]
	-m,	mode of editing, range 0-2, [default=0]
			0: best substitution, or first good indel
			1: best substitution, or best indel
			2: best edit overall (suggestion that you reduce i and d for performance)
	-v,	verbose mode (-v 1 = yes, default = 0, no)

	--help,		display this message and exit 
	--version,	output version information and exit

</pre>

### ntEdit Modes
------------
<pre>
Mode 0: (default)
	ntEdit will try to substitute the last base of an incorrect k-mer with a different ATGC base. If that k-mer is found in the bloom filter and has enough subset support, ntEdit will then try the other substitution bases and then choose the best substitution fix. However, if the substituion was not found, then ntEdit will try all indels of max length (-i) and (-d) starting with that substitution base and make edit based on the first accepted indel. 

Mode 1: 
	ntEdit will try to substitute the last base of an incorrect k-mer with a different ATGC base. If that k-mer is found in the bloom filter and has enough subset support, ntEdit will then try the other substitution bases and then choose the best substitution fix. However, if the substitution was nto found,t hen ntEdit will try all indels of max length (-i) and (-d) starting witht hat substitution base and make edit based on the best accepted indel. 

Mode 2: 
	ntEdit will choose the best substitution or indel for each incorrect k-mer. Since this can be very computationally expensive because ntEdit tries every combination possible, it is recommended that you reduce (-i) and (-d). 
</pre>

*We recommend running ntEdit in Mode 1 (or 0)


### Test data
---------
<pre>
Go to ./demo
</pre>
(cd demo)


run:
-------------------------------------
./runme.sh

(../ntedit -f ecoliWithMismatches001Indels0001.fa.gz -k 25 -r solidBF_k25.bf -d 5 -i 4 -b ntEditEcolik25)

ntEdit will polish an E. coli genome sequence with substitution error ~0.001 and indels ~0.0001 using pre-made ntHits Bloom filter

Expected files will be:
ntEditEcolik25_changes.tsv
ntEditEcolik25_edited.fa

Compare with:
nteditk25_changes.tsv
nteditk25_edited.fa


### How it works
------------
![Logo](https://github.com/bcgsc/ntEdit/blob/master/figS1.png)
Sequence reads are first shredded into kmers using ntHits, keeping track of kmer multiplicity. The kmers that pass coverage thresholds (ntHits, -c option builds a filter with kmers having a coverage higher than c) are used to construct a Bloom filter (BF). The draft assembly is supplied to ntEdit (-f option, fasta file), along with the BF (-r option) and sequences are read sequentially. Sequence strings are shredded into words of length k (kmers) at a specified value (-k option) matching that used to build the BF, and each kmer from 5’ to 3’ queries the BF data structure for presence/absence (step 1). When a kmer is not found in the filter, a subset (Sk) of overlapping k kmers (defined by k over three, k/3) containing the 3’-end base is queried for absence (step 2). The subset Sk, representing a subsampling of k kmers obtained by sliding over 3 bases at a time over k bases, is chosen to minimize the number of checks against the Bloom filter. Of this subset, when the number of absent (-) kmers matches or exceeds a threshold defined by Sk- >= k/x (-x option), representing the majority of kmers in Sk, editing takes place (step 3 and beyond), otherwise step 1 resumes. In the former case, the 3’-end base is permuted to one of the three alternate bases (step 3), and the subset (Sk_alt) containing the change is assessed for Bloom filter presence (+). When that number matches or exceeds the threshold defined by Sk_alt+ >= k/y (-y option), which means the base substitution qualifies, it is tracked along with the number of supported kmers and the remaining alternate 3’-end base substitutions are also assessed (ie. resuming step 3 until all bases inspected). If the edit does not qualify, then a cycle of base insertion(s) and deletion(s) of up to –i and –d bases begins (step 4, -i option and step 5, -d option, respectively). As is the case for the substitutions, a subset of k kmers containing the indel change is tested for presence. If there are no qualifying changes, then the next alternate 3’-end base is inspected as per above; otherwise the change is applied to the sequence string and the next assembly kmer is inspected (step 1). The process is repeated until a qualifying change or until no suitable edits are found. In the latter case, we go back to step 1. When a change is made, the position on the new sequence is tracked, along with an alternate base with lesser or equal k kmer subset support, when applicable. Currently, ntEdit only tracks cases when edits are made (steps 3-5), and does not flag unedited, missing draft kmers (steps 1-2).  


### OUTPUT FILES
------------------------

|Output files|                    Description|
|---|---|
|_changes.tsv                 | tab-separated file; ID      bpPosition+1    OriginalBase    NewBase Support 25-mers (out of k/3)   AlternateNewBase   Alt.Support k-mers   eg. U00096.3_MG1655_k12     117     A       T       9|
|_edited.fa                   | fasta file; contains the polished genome assembly |



### License
-------

ntEdit Copyright (c) 2018-2019 British Columbia Cancer Agency Branch.  All rights reserved.

ntEdit is released under the GNU General Public License v3

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, version 3.
 
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.
 
You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.
 
For commercial licensing options, please contact
Patrick Rebstein <prebstein@bccancer.bc.ca>
