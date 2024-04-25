#!/usr/bin/env bash

set -e

echo DOWNLOADING DATA
wget -nc https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_1.fq.gz https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_2.fq.gz

run-ntedit polish -k 25 -d 5 -i 4 --reads Sim_100_300 --draft ecoliWithMismatches001Indels0001.fa.gz
diff -q ntedit_k25_edited.fa ecoli_ntedit_k25_edited.fa 
diff -q ntedit_k25_changes.tsv ecoli_ntedit_k25_changes.tsv
