#!/usr/bin/env snakemake -s

# Snakefile for ntedit pipeline
import os

# Read parameters from config or set default values
draft=config["draft"]
reads_prefix=config["reads"] if "reads" in config else ""
k=config["k"]

# ancestry parameters
genomes = config["genomes"] if "genomes" in config else ""
genome_prefix = ".".join([os.path.basename(os.path.realpath(genome)).removesuffix(".fa").removesuffix(".fasta").removesuffix(".fna") for genome in genomes])

if draft is None or reads_prefix is None or k is None:
    raise ValueError("You must specify draft, reads, and k in your command. See 'snakemake -s ntedit_run_pipeline.smk help' for more information.")
# Check that draft file exists
if not os.path.isfile(draft):
    raise ValueError("Draft file does not exist. Please check that the file name is correct and that the file is in the current working directory.")
# Check that reads files exist
if [file for file in os.listdir('.') if file.startswith(reads_prefix) and file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))] == []:
    if genomes == "":
        raise ValueError("Reads files do not exist. Please check that the prefix is correct and that the files are in the current working directory.")

reads_files = [file for file in os.listdir('.') if file.startswith(reads_prefix) and file.endswith((".fq", ".fq.gz", ".fastq", ".fastq.gz"))]

# Check that k is an integer
try:
    k = int(k)
except ValueError:
    raise ValueError("k must be an integer.")
# Check that k is greater than 0
if k <= 0:
    raise ValueError("k must be greater than 0.")

# Common parameters
t = config["t"] if "t" in config else 4
b = config["b"] + "_" if "b" in config and config["b"] != "" else ""

# ntHits parameters
solid = config["solid"] if "solid" in config else False
cutoff = config["cutoff"] if "cutoff" in config else 2

# ntEdit parameters
z = config["z"] if "z" in config else 100
i = config["i"] if "i" in config else 5
d = config["d"] if "d" in config else 5
x = config["x"] if "x" in config else 5.000
y = config["y"] if "y" in config else 9.000
cap = config["cap"] if "cap" in config else (int(k * 3 / 2) if k is not None else None)
m = config["m"] if "m" in config else 0
v = config["v"] if "v" in config else 0
a = config["a"] if "a" in config else 0
j = config["j"] if "j" in config else 3
s = config["s"] if "s" in config else 0
X = config["X"] if "X" in config else -1
Y = config["Y"] if "Y" in config else -1
p = config["p"] if "p" in config else 1
q = config["q"] if "q" in config else 255
l = config["l"] if "l" in config else ""


if l != "":
    if not os.path.isfile(l):
        raise ValueError("VCF file not found. Please check that the file name is correct.")

if X != -1 or Y != -1:
    if X == -1:
        X = 0.5
    if Y == -1:
        Y = 0.5

# time command
mac_time_command = "command time -l -o"
linux_time_command = "command time -v -o"
time_command = mac_time_command if os.uname().sysname == "Darwin" else linux_time_command

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
        echo "    k        k-mer size, REQUIRED"
        echo "    t        number of threads [default=4]"
        echo "    b        output file prefix, OPTIONAL"
        echo ""
        echo "Options specific to ntHits:"
        echo "    solid    output the solid k-mers (non-erroneous k-mers), True = yes, False = no [default=False]"
        echo "	  cutoff   the minimum coverage of k-mers in output bloom filter, [default=2, ignored if solid=True]"
        echo ""
        echo "Options specific to ntEdit:"
        echo "    z        minimum contig length [default=100]"
        echo "    i        maximum number of insertion bases to try, range 0-5, [default=5]"
        echo "    d        maximum number of deletions bases to try, range 0-10, [default=5]"
        echo "    x        k/x ratio for the number of k-mers that should be missing, [default=5.000]"
        echo "    y        k/y ratio for the number of edited k-mers that should be present, [default=9.000]"
        echo "    j        controls size of k-mer subset. When checking subset of k-mers, check every jth k-mer, [default=3]"
        echo "    cap      cap for the number of base insertions that can be made at one position, [default=k*1.5]"
        echo "    X        ratio of number of k-mers in the k subset that should be missing in order to attempt fix (higher=stringent), [default=0.5]"
        echo "    Y        ratio of number of k-mers in the k subset that should be present to accept an edit (higher=stringent), [default=0.5]"
        echo "    m        mode of editing, range 0-2, [default=0]"
        echo "                 0: best substitution, or first good indel"
        echo "                 1: best substitution, or best indel"
        echo "                 2: best edit overall (suggestion that you reduce i and d for performance)"
        echo "    a        Soft masks missing k-mer positions having no fix (-v 1 = yes, default = 0, no)"
        echo "    s        SNV mode. Overrides draft k-mer checks, forcing reassessment at each position (1 = yes, default = 0, no. EXPERIMENTAL)"
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
        edited_draft=f"{b}ntedit_k{k}_edited.fa"
    params:
        prefix = f"{b}ntedit_k{k}",
        benchmark = f"{time_command} ntedit_{reads_prefix}.time",
        ratio = f"-X {X} -Y {Y}" if X != -1 or Y != -1 else "",
        vcf = f"-l {l}" if l != "" else ""
    shell:
        "{params.benchmark} ntedit -r {input.bloom_filter} -f {draft} -b {params.prefix} -t {t} -z {z} -i {i} -d {d} -x {x} -y {y} -c {cap} -m {m} -v {v} -a {a} -j {j} {params.ratio} -s {s} {params.vcf}"
        

rule nthits:
    input:
        hist=f"{reads_prefix}_k{k}.hist",
        reads_files = reads_files
    output:
        bloom_filter=f"{reads_prefix}_k{k}.bf"
    params:
        benchmark = f"{time_command} nthits_{reads_prefix}_k{k}.time",
        min_cutoff = f"--solid" if solid else f"-cmin {cutoff}"      
    shell:
        "{params.benchmark} nthits bf -t {t} -f {input.hist} -o {output.bloom_filter} -k {k} {params.min_cutoff} -v {input.reads_files}"
        

rule ntcard:
    input:
        reads_files = reads_files
    output:
        hist=f"{reads_prefix}_k{k}.hist"
    params:
        benchmark = f"{time_command} ntcard_{reads_prefix}_k{k}.time"
    shell:
        "{params.benchmark} ntcard -k {k} -t {t} -p {reads_prefix} {input.reads_files}"


rule ntedit_cbf:
    input:
        bloom_filter=f"{reads_prefix}_k{k}.cbf"
    output:
        edited_draft=f"{b}ntedit_k{k}_cbf_p{p}_q{q}_edited.fa"
    params:
        prefix = f"{b}ntedit_k{k}_cbf_p{p}_q{q}",
        benchmark = f"{time_command} ntedit_{reads_prefix}_k{k}_cbf_p{p}_q{q}.time",
        ratio = f"-X {X} -Y {Y}" if X != -1 or Y != -1 else "",
        vcf = f"-l {l}" if l != "" else ""
    shell:
        "{params.benchmark} ntedit -r {input.bloom_filter} -f {draft} -b  {params.prefix} -p {p} -q {q} -t {t} -z {z} -i {i} -d {d} -x {x} -y {y} -c {cap} -m {m} -v {v} -a {a} -j {j} {params.ratio} -s {s} {params.vcf}"        

rule nthits_cbf:
    input:
        hist=f"{reads_prefix}_k{k}.hist",
        reads_files = reads_files
    output:
        bloom_filter=f"{reads_prefix}_k{k}.cbf"
    params:
        benchmark = f"{time_command} nthits_{reads_prefix}_k{k}_cbf.time",
        min_cutoff = f"--solid" if solid else f"-cmin {cutoff}"      
    shell:
        "{params.benchmark} nthits cbf -t {t} -f {input.hist} -o {output.bloom_filter} -k {k} {params.min_cutoff} {input.reads_files}"


rule ntedit_snv_genome:
    input: f"{genome_prefix}_ntedit_k{k}_variants.vcf"

rule ntedit_snv_reads:
    input: f"{reads_prefix}_ntedit_k{k}_variants.vcf"

rule ntedit_snv_reads_cbf:
    input: f"{reads_prefix}_ntedit_k{k}_cbf_variants.vcf"

rule ntedit_snv_cbf:
    input:
        bloom_filter = expand("{{prefix}}_k{k}.cbf", k=k)
    output:
        out_vcf=expand("{{prefix}}_ntedit_k{k}_cbf_variants.vcf", k=k)
    params:
        prefix = expand("{{prefix}}_ntedit_k{k}_cbf", k=k),
        benchmark = expand("{time_command} ntedit_{{prefix}}.time", time_command=time_command),
        ratio = f"-X {X} -Y {Y}" if X != -1 or Y != -1 else "",
        vcf = f"-l {l}" if l != "" else ""
    shell:
        "{params.benchmark} ntedit -r {input.bloom_filter} -f {draft} -b {params.prefix} -t {t} -z {z} -y {y} -v {v} -a {a} -j {j} {params.ratio} -s 1 {params.vcf}"

rule ntedit_snv:
    input:
        bloom_filter = expand("{{prefix}}_k{k}.bf", k=k)
    output:
        out_vcf=expand("{{prefix}}_ntedit_k{k}_variants.vcf", k=k)
    params:
        prefix = expand("{{prefix}}_ntedit_k{k}", k=k),
        benchmark = expand("{time_command} ntedit_{{prefix}}.time", time_command=time_command),
        ratio = f"-X {X} -Y {Y}" if X != -1 or Y != -1 else "",
        vcf = f"-l {l}" if l != "" else ""
    shell:
        "{params.benchmark} ntedit -r {input.bloom_filter} -f {draft} -b {params.prefix} -t {t} -z {z} -y {y} -v {v} -a {a} -j {j} {params.ratio} -s 1 {params.vcf}"


rule ntedit_genome_bf:
    input:
        ntcard_f_stats = f"{genome_prefix}.k{k}.tsv",
        genomes = genomes
    output:
        f"{genome_prefix}_k{k}.bf"
    params:
        options = f"-k {k} -t {t}",
        benchmark = f"{time_command} make_bf_{genome_prefix}_k{k}.time"
    run:
        with open(input.ntcard_f_stats, 'r', encoding="utf8") as fin:
            for line in fin:
                _, f, num = line.strip().split("\t")
                if f == "F0":
                    num_elements = num
                    break
        shell("{params.benchmark} make_genome_bf --genome {input.genomes} {params.options} --num_elements {num_elements} -o {output}")
        

rule genomes_ntcard:
    input: 
        genomes = genomes
    output:
        hist=f"{genome_prefix}.k{k}.hist",
        f_stats=f"{genome_prefix}.k{k}.tsv"
    params:
        benchmark = f"{time_command} ntcard_{genome_prefix}_k{k}.time",
        options = f"-t {t} -k {k}"
    shell:
        "{params.benchmark} ntcard {params.options} -o {output.hist} {input.genomes} &> {output.f_stats}"
