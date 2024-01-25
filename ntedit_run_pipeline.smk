#!/usr/bin/env snakemake -s

# Snakefile for ntedit pipeline
import os

# Read parameters from config or set default values
draft=config["draft"]
reads_prefix=config["reads"]
k=config["k"]

if draft is None or reads_prefix is None or k is None:
    raise ValueError("You must specify draft, reads, and k in your command. See 'snakemake -s ntedit_run_pipeline.smk help' for more information.")
# Check that draft file exists
if not os.path.isfile(draft):
    raise ValueError("Draft file does not exist. Please check that the file name is correct and that the file is in the current working directory.")
# Check that reads files exist
if [file for file in os.listdir('.') if file.startswith(config["reads"]) and file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))] == []:
    raise ValueError("Reads files do not exist. Please check that the prefix is correct and that the files are in the current working directory.")
# Check that k is an integer
try:
    k = int(k)
except ValueError:
    raise ValueError("k must be an integer.")
# Check that k is greater than 0
if k <= 0:
    raise ValueError("k must be greater than 0.")

# Common parameters
t = config["t"] if "t" in config else 1
b = config["b"] + "_" if "b" in config and config["b"] != "" else ""

# ntHits parameters
solid = config["solid"] if "solid" in config else False
cutoff = config["cutoff"] if "cutoff" in config else 2

# ntEdit parameters
z = config["z"] if "z" in config else 100
i = config["i"] if "i" in config else 4
d = config["d"] if "d" in config else 5
x = config["x"] if "x" in config else 5.000
y = config["y"] if "y" in config else 9.000
cap = config["cap"] if "cap" in config else (int(k * 3 / 2) if k is not None else None)
m = config["m"] if "m" in config else 0
v = config["v"] if "v" in config else 0
a = config["a"] if "a" in config else 0
j = config["j"] if "j" in config else 3
s = config["s"] if "s" in config else 0
X = config["X"] if "X" in config else 0.5
Y = config["Y"] if "Y" in config else 0.5
p = config["p"] if "p" in config else 1
q = config["q"] if "q" in config else 255

# time command
time_command = "command time -v -o"

rule all:
    input:
        f"{b}ntedit_k{k}_edited.fa"

rule help:
    shell:
        """
        echo ""
        echo "Usage: snakemake -s ntedit_run_pipeline.smk [target]"
        echo ""
        echo "Targets:"
        echo "    ntedit   run ntedit to polish draft genome assembly, REQUIRED"
        echo "    ntedit_cbf   run ntedit using a counting Bloom filter to polish draft genome assembly, REQUIRED"
        echo ""
        echo "Options:"
        echo "    draft    draft genome assembly. Must be specified with exact FILE NAME. Ex: draft=myDraft.fa (FASTA, Multi-FASTA, and/or gzipped compatible), REQUIRED"
        echo "    reads    prefix of reads file(s). All files in the working directory with the specified prefix will be used for polishing (fastq, fasta, gz, bz, zip), REQUIRED"
        echo "    time     logs time and memory usage to file for main steps (Set to 1 to enable logging)"
        echo "    k        kmer size, REQUIRED"
        echo "    t        number of threads [default=1]"
        echo "    b        output file prefix, OPTIONAL"
        echo ""
        echo "Options specific to ntHits:"
        echo "    solid    output the solid k-mers (non-erroneous k-mers), True = yes, False = no [default=False]"
        echo "	  cutoff   the minimum coverage of kmers in output bloom filter, [default=2, ignored if solid=True]"
        echo ""
        echo "Options specific to ntEdit:"
        echo "    z        minimum contig length [default=100]"
        echo "    i        maximum number of insertion bases to try, range 0-5, [default=4]"
        echo "    d        maximum number of deletions bases to try, range 0-5, [default=5]"
        echo "    x        k/x ratio for the number of kmers that should be missing, [default=5.000]"
        echo "    y        k/y ratio for the number of edited kmers that should be present, [default=9.000]"
        echo "    j        controls size of kmer subset. When checking subset of kmers, check every jth kmer, [default=3]"
        echo "    cap      cap for the number of base insertions that can be made at one position, [default=k*1.5]"
        echo "    X        ratio of number of kmers in the k subset that should be missing in order to attempt fix (higher=stringent), [default=0.5]"
        echo "    Y        ratio of number of kmers in the k subset that should be present to accept an edit (higher=stringent), [default=0.5]"
        echo "    m        mode of editing, range 0-2, [default=0]"
        echo "                 0: best substitution, or first good indel"
        echo "                 1: best substitution, or best indel"
        echo "                 2: best edit overall (suggestion that you reduce i and d for performance)"
        echo "    a        Soft masks missing kmer positions having no fix (-v 1 = yes, default = 0, no)"
        echo "    s        SNV mode. Overrides draft kmer checks, forcing reassessment at each position (1 = yes, default = 0, no. EXPERIMENTAL)"
        echo "    v        verbose mode (1 = yes, default = 0, no)"
        echo "    p, minimum k-mer coverage threshold (CBF only) [default=minimum of counting Bloom filter\n"
        echo "    counts, cannot be larger than 255]\n"
        echo "    q, maximum k-mer coverage threshold (CBF only) [default=255, largest possible value]\n"
        echo ""
        echo "Example: Polishing myDraft.fa with myReads1.fq and myReads2.fq"
        echo "        snakemake -s ntedit_run_pipeline.smk ntedit draft=myDraft.fa reads=myReads cutoff=2 or"
        echo "        snakemake -s ntedit_run_pipeline.smk ntedit draft=myDraft.fa reads=myReads solid=true"
        echo ""
        echo "Make sure your read files all have the same prefix, as indicated by 'reads=<prefix>'. The Snakefile will use all files in the current working directory with this prefix for polishing."
        echo "To ensure that the pipeline runs correctly, make sure that the following tools are in your PATH: ntedit, nthits"
	    echo "You must either specify the cutoff parameter or define solid=true in your command to set it automatically"
	    echo "If one of X/Y is set, ntEdit will use those parameters instead. Otherwise, it uses x/y by default."
        """

        
rule ntedit:
    input:
        bloom_filter=f"{reads_prefix}_k{k}.bf"
    output:
        edited_draft=f"ntedit_k{k}_edited.fa"
    run:
        shell(
            f"{time_command} ntedit_{reads_prefix}_k{k}.time ntedit -r {input.bloom_filter} -f {input.draft} -b {b}ntedit_k{k} -t {t} -z {z} -i {i} -d {d} -x {x} -y {y} -c {cap} -m {m} -v {v} -a {a} -j {j} -X {X} -Y {Y} -s {s}"
        )

rule nthits:
    input:
        hist=f"{reads_prefix}_k{k}.hist",
        reads_files = [file for file in os.listdir('.') if file.startswith(config["reads"]) and file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))]
    output:
        bloom_filter=f"{reads_prefix}_k{k}.bf"        
    run:
        shell(
            f"{time_command} nthits_{reads_prefix}_k{k}.time nthits bf -t {t} -f {input.hist} -o {output.bloom_filter} -k {k} {('--solid ' if solid else '')}-cmin {cutoff} {input.reads_files}"
        )

rule ntcard:
    input:
        reads_files = [file for file in os.listdir('.') if file.startswith(config["reads"]) and file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))]
    output:
        hist=f"{reads_prefix}_k{k}.hist"
    run:
        shell(
            f"{time_command} ntcard_{reads_prefix}_k{k}.time ntcard -k {k} -t {t} -o {output.hist} {input.reads_files}"
        )

rule ntedit_cbf:
    input:
        bloom_filter=f"{reads_prefix}_k{k}.cbf"
    output:
        edited_draft=f"ntedit_k{k}_cbf_p{p}_q{q}_edited.fa"
    run:
        shell(
            f"{time_command} ntedit_{reads_prefix}_k{k}_cbf_p{p}_q{q}.time ntedit -r {input.bloom_filter} -f {input.draft} -b {b}ntedit_k{k}_cbf_p{p}_q{q} -p {p} -q {q} -t {t} -z {z} -i {i} -d {d} -x {x} -y {y} -c {cap} -m {m} -v {v} -a {a} -j {j} -X {X} -Y {Y} -s {s}"
        )

rule nthits_cbf:
    input:
        hist=f"{reads_prefix}_k{k}.hist",
        reads_files = [file for file in os.listdir('.') if file.startswith(config["reads"]) and file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))]
    output:
        bloom_filter=f"{reads_prefix}_k{k}.cbf"        
    run:
        shell(
            f"{time_command} nthits_{reads_prefix}_k{k}_cbf.time nthits cbf -t {t} -f {input.hist} -o {output.bloom_filter} -k {k} {('--solid ' if solid else '')}-cmin {cutoff} {input.reads_files}"
        )
