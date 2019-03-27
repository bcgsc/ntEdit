![Logo](https://github.com/bcgsc/ntEdit/blob/master/ntedit-logo.png)

# ntEdit

## Scalable genome assembly polishing
## 12/2018
## email: rwarren [at] bcgsc [dot] ca


### Description
-----------

ntEdit is a genomics application for polishing genome assembly drafts.
ntEdit simplifies polishing and "haploidization" of gene and genome sequences with its re-usable Bloom filter design.
We expect ntEdit to have additional applications in fast mapping of simple nucleotide variations between any two individuals or species’ genomes.


### Implementation and requirements
-------------------------------

ntEdit v1.0 is prototyped in PERL and runs on any OS where PERL is installed.


### Install
-------

Download the tar ball, gunzip and extract the files on your system using:
<pre>
gunzip ntedit_v1-0-1.tar.gz
tar -xvf ntedit_v1-0-1.tar
</pre>

### Dependencies
-------

1. ntHits (https://github.com/bcgsc/nthits)
2. BloomFilter utilities (provided in ./lib - you may need to recompile, see instructions below)


### Instructions for building the nthits Bloom filter utility PERL module
-------

1. BUILD a PERL5 module

Make sure you have swig installed and included in your path.

http://www.swig.org/


TO BUILD a Perl5 module (run in swig/):
```
a) /home/<user>/<path to swig>/preinst-swig -Wall -c++ -perl5 BloomFilter.i

b) g++ -c BloomFilter_wrap.cxx -I/usr/lib64/perl5/CORE -fPIC -Dbool=char -O3

The location of CORE may change on your system.
If you use linuxbrew, this works:
-I/<yourpath>/linuxbrew/Cellar/perl/5.28.0/lib/perl5/5.28.0/x86_64-linux-thread-multi/CORE/

c) g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so -O3
```

TO COMPILE, swig needs the following Perl5 headers:
```C++
#include "EXTERN.h"
#include "perl.h"
#include "XSUB.h"
```
If they are not located in /usr/lib64/perl5, you can run "perl -e 'use Config; print $Config{archlib};" to locate them.


2. VERIFY your install

in the lib folder, execute:
./test.pl

-All tests should pass


3. CHANGE the relative path to BloomFilter.pm in ntEdit.pl/test.pl 

You only need to change if you have re-built in a relative directory different
from:
<pre>
use lib "$FindBin::Bin/./lib/"; (for LINKS and test.pl)
</pre>


### Documentation
-------------

Refer to the README.md file on how to run ntEdit and our manuscript for information about the software and its performance 
Questions or comments?  We would love to hear from you!
rwarren at bcgsc.ca


### Citing ntEdit 
------------

<pre>
René L Warren, Lauren Coombe, Hamid Mohamadi, Jessica Zhang, Barry Jaquish, Nathalie Isabel, Steven JM Jones, Jean Bousquet, Joerg Bohlmann and Inanç Birol
ntEdit: scalable genome assembly polishing
TBD
</pre>

The experimental data described in our paper can be downloaded here: http://www.bcgsc.ca/downloads/btl/ntedit/
Thank you for using, developing and promoting this free software.


### Credits
-------

ntEdit:
Rene Warren

nthits / nthash / BloomFilter.pm:
Hamid Mohamadi (recursive/ntHash)

C++ implementation
Jessica Zhang


### How to run in a pipeline
-------

1. Running nthits (please see nthits documentation)
nthits -c <kmer coverage threshold> -k <kmer length> -j <number of threads> reads
eg.
./nthits -c 1 --outbloom --solid -k 25 -j 48 Sim_100_300_1.fq Sim_100_300_2.fq
or
./nthits -c 1 --outbloom --solid -k 25 -j 48 @reads.in

Where @reads.in is a file listing the path to all fastq files


2. Running ntEdit (see complete usage below)
ntEdit.pl -f <fasta file to polish> -k <kmer length> -r <Bloom filter from nthits>
eg.
./ntEdit.pl -f ecoliWithMismatches001Indels0001.fa -r solidBF_k25.bf -k 25 -b ntEditEcolik25


### Running ntEdit
-------------
<pre>
e.g. ./ntEdit.pl -f ecoliWithMismatches001Indels0001.fa -r solidBF_k25.bf -k 25 -b ntEditEcolik25

Usage: ../ntEdit.pl [v1.0.1]
-f  draft genome assembly (Multi-FASTA format, required)
-r  Bloom filter of sequence reads (ntHits format, required)
-k  k-mer value (default -k 35, optional, same value of k to build -r Bloom filter)
-d  maximum number of base deletions (default -d 0, optional, range 1-5)
-i  maximum number of base insertions (default -i 0, optional, range 1-5, values higher than 4 will impact run speed)
-x  leniency factor 1 : determines whether a missing kmer should be edited (default -x 5, optional, lower value=less permissive)
-y  leniency factor 2 : determines whether a change should be kept (default -y 9, optional, lower value=less permissive)
-z  minimum contig length to consider (default -z 100, optional)
-b  base name for your output files (optional)
-v  Runs in verbose mode (-v 1 = yes, default = no, optional)

NOTE: BLOOM FILTER MUST BE DERIVED FROM THE SAME SEQUENCE READS FILE(S) SUPPLIED to ntHits, WITH SAME -k VALUE SPECIFIED

Error: Missing mandatory options -f, -r

</pre>


### Test data
---------
<pre>
Go to ./demo
</pre>
(cd demo)


run:
-------------------------------------
./runme.sh (../ntEdit.pl -f ecoliWithMismatches001Indels0001.fa.gz -r solidBF_k25.bf -k 25 -b ntEditEcolik25)

ntEdit will polish an E. coli genome sequence with substitution error ~0.001 and indels ~0.0001 using pre-made nthits Bloom filter

Expected files will be:
ntEditEcolik25.log
ntEditEcolik25_changes.tsv
ntEditEcolik25_edited.fa


### How it works
------------
![Logo](https://github.com/bcgsc/ntEdit/blob/master/figS1.png)
Sequence reads are first shredded into kmers using ntHits, keeping track of kmer multiplicity. The kmers that pass coverage thresholds (ntHits, -c option) are used to construct a Bloom filter (BF). The draft assembly is supplied to ntEdit (-f option, fasta file), along with the BF (-r option) and sequences are read sequentially. Draft assembly contigs are shredded into kmers (at a specified –k value matching that used to build the BF), and each kmer from 5’ to 3’ queries the BF data structure for presence/absence (step 1). When kmers are not found in the filter, a subset (k/3) of k kmers containing the 3’-end base is queried for absence (step 2). When >=k/5 kmers (by default, -x option) are absent, editing takes place, otherwise step 1 resumes and the next assembly kmer is assessed. The 3’-end base is permuted to one of the three alternate bases (step 3), and a subset of k kmers containing the change is assessed (>= k/9 by default, -y option). When a base substitution is made that qualifies, it is tracked along with the number of supported kmers and the remaining alternate 3’-end base substitutions are also assessed (ie. resuming step3 until all bases inspected). If it does not qualify, then a cycle of base insertion(s) (step 4) and deletion(s) (step 5) begins. As is the case for the substitutions, a subset of k kmers containing the indel change is assessed (>= k/9 by default, -y option). If there are no qualifying changes, then the next alternate 3’-end base is inspected as per above; otherwise the change is applied to the sequence string and the next assembly kmer is inspected (step 1). The process is repeated until a qualifying change or until no suitable edits are found (back to step 1).  

### OUTPUT FILES
------------------------

|Output files|                    Description|
|---|---|
|.log                         | text file; Logs execution time / errors / pairing stats|
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
