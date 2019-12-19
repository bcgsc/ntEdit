#source /home/hmohamadi/.e7_env
#RLW 03/2019
#echo Running ntHits
#==================
#wget https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_1.fq.gz
#wget https://www.bcgsc.ca/downloads/btl/ntedit/sim_reads_ecoli/Sim_100_300_2.fq.gz
#/usr/bin/time -v -o nthitk25.time nthits -c 2 -b 36 -k 25 --outbloom -p solidBF Sim_100_300_1.fq.gz Sim_100_300_2.fq.gz
echo Running ntEdit...
#===================
/usr/bin/time -v -o nteditk25.time ../ntedit -f ecoliWithMismatches001Indels0001.fa.gz -r solidBF_k25.bf -d 5 -i 4 -b nteditk25
echo Output corrections in nteditk25_
