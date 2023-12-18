#!/usr/bin/env bash

set -e

echo DOWNLOADING DATA
wget -nc https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_1.fq.gz https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_2.fq.gz

echo -e "\nRUNNING NTCARD"
ntcard -k 25 -c 64 -p ntcard Sim_100_300_1.fq.gz Sim_100_300_2.fq.gz

echo -e "\nRUNNING NTHITS"
nthits cbf -v -f ntcard_k25.hist -o nthits.cbf -k 25 --solid Sim_100_300_1.fq.gz Sim_100_300_2.fq.gz
# /usr/bin/time -v -o nthits.time nthits bf -v -f ntcard_k25.hist -o nthits.bf -k 25 --solid Sim_100_300_1.fq.gz Sim_100_300_2.fq.gz

echo -e "\nRUNNING NTEDIT"
ntedit -f ecoliWithMismatches001Indels0001.fa.gz -r nthits.cbf -b ntedit
