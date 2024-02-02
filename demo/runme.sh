#!/usr/bin/env bash

set -e

echo DOWNLOADING DATA
wget -nc https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_1.fq.gz https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_2.fq.gz

run-ntedit --bloomType cbf -k 31 --reads Sim_100_300 --draft ecoliWithMismatches001Indels0001.fa.gz
