#!/usr/bin/env python3
'''
ntEdit: Fast, lightweight, scalable genome sequence polishing and SNV detection & annotation
Written by Johnathan Wong, adapted from https://github.com/bcgsc/ntSynt/blob/main/bin/ntSynt
'''
import argparse
import os
import shlex
import subprocess
from packaging import version

import snakemake

NTEDIT_VERSION = "ntEdit v2.0.0"

def main():
    "Run ntEdit snakemake file"

    parser = argparse.ArgumentParser(description="ntEdit: Fast, lightweight, scalable genome sequence polishing and SNV detection & annotation",
                                     formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument("--mode",
                        help="Mode of ntEdit (bf, cbf), [default=bf]", default="bf", type=str)    

    parser.add_argument("-k",
                        help="Kmer size, REQUIRED",
                        required=True, type=int)

    parser.add_argument("--draft",
                        help="Draft genome assembly. Must be specified with exact FILE NAME. "
                            "Ex: --draft myDraft.fa (FASTA, Multi-FASTA, and/or gzipped compatible), REQUIRED",
                        required=True)

    parser.add_argument("--reads",
                        help="Prefix of reads file(s). All files in the working directory with the specified prefix will be used for polishing "
                            "(fastq, fasta, gz, bz, zip), REQUIRED",
                        required=True)

    parser.add_argument("--cutoff",
                        help="The minimum coverage of kmers in output Bloom filter, [default=2, ignored if solid=True]",
                        default=2, type=int)

    parser.add_argument("-t",
                        help="Number of threads [default=1]", default=1, type=int)

    parser.add_argument("--solid",
                        help="Output the solid k-mers (non-erroneous k-mers), True = yes, False = no [default=False]",
                        action="store_true", default=False)

    parser.add_argument("-z",
                        help="Minimum contig length [default=100]", default=100, type=int)

    parser.add_argument("-i",
                        help="Maximum number of insertion bases to try, range 0-5, [default=4]", default=4, type=int)

    parser.add_argument("-d",
                        help="Maximum number of deletions bases to try, range 0-5, [default=5]", default=5, type=int)

    parser.add_argument("-x",
                        help="K/x ratio for the number of kmers that should be missing, [default=5.000]", default=5.000, type=float)

    parser.add_argument("-y",
                        help="K/y ratio for the number of edited kmers that should be present, [default=9.000]", default=9.000, type=float)

    parser.add_argument("--cap",
                        help="Cap for the number of base insertions that can be made at one position, [default=k*1.5]",
                        default=None, type=float)

    parser.add_argument("-m",
                        help="Mode of editing, range 0-2, [default=0]", default=0, type=int)

    parser.add_argument("-v",
                        help="Verbose mode (1 = yes, default = 0, no)", default=0, type=int)

    parser.add_argument("-X",
                        help="Ratio of number of kmers in the k subset that should be missing in order to attempt fix (higher=stringent), "
                            "[default=0.5]", default=0.5, type=float)

    parser.add_argument("-Y",
                        help="Ratio of number of kmers in the k subset that should be present to accept an edit (higher=stringent), "
                            "[default=0.5]", default=0.5, type=float)

    parser.add_argument("-p",
                        help="Minimum k-mer coverage threshold (CBF only) [default=1]", default=1, type=int)

    parser.add_argument("-q",
                        help="Maximum k-mer coverage threshold (CBF only) [default=255, largest possible value]", default=255, type=int)
    
    parser.add_argument("-V", "--version", action="version", version=NTEDIT_VERSION)

    parser.add_argument("-n", "--dry-run", help="Print out the commands that will be executed", action="store_true")

    parser.add_argument("-f", "--force", help="Run all ntEdit steps, regardless of existing output files",
                        action="store_true")

    args = parser.parse_args()

    base_dir = os.path.dirname(os.path.realpath(__file__))

    if args.mode == "bf":
        mode = "ntedit"
    elif args.mode == "cbf":
        mode = "ntedit_cbf"
    else:
        raise ValueError(f"Invalid mode {args.mode}")
    

    intro_string = ["Running ntedit...",
                    f"Parameter settings:",
                    f"\t--mode {mode}",
                    f"\t--k {args.k}",
                    f"\t--draft {args.draft}",
                    f"\t--reads {args.reads}",
                    f"\t--cutoff {args.cutoff}",
                    f"\t--t {args.t}",
                    f"\t--solid {args.solid}",
                    f"\t--z {args.z}",
                    f"\t--i {args.i}",
                    f"\t--d {args.d}",
                    f"\t--x {args.x}",
                    f"\t--y {args.y}",
                    f"\t--cap {args.cap}",
                    f"\t--m {args.m}",
                    f"\t--v {args.v}",
                    f"\t--X {args.X}",
                    f"\t--Y {args.Y}",
                    f"\t--p {args.p}",
                    f"\t--q {args.q}",
                    ]

    print("\n".join(intro_string), flush=True)

    command = f"snakemake -s {base_dir}/ntedit_run_pipeline.smk {mode} -p --cores {args.t} " \
            f"--config draft={args.draft} reads={args.reads} k={args.k} cutoff={args.cutoff} t={args.t} " \
            f"solid={args.solid} z={args.z} i={args.i} d={args.d} x={args.x} y={args.y} " \
            f"m={args.m} v={args.v} X={args.X} Y={args.Y} p={args.p} q={args.q} "

    if args.cap is not None:
        command += "cap={args.cap}"

    if version.parse(snakemake.__version__) >= version.parse("7.8.0"): # Keep behaviour consistent for smk versions
        command += "--rerun-trigger mtime "

    if args.dry_run:
        command += " -n"

    if args.force:
        command += " -F"

    print(f"Running {command}", flush=True)

    command = shlex.split(command)

    ret = subprocess.call(command)
    if ret != 0:
        raise subprocess.SubprocessError("ntEdit failed - check the logs for the error.")

    print("Done ntEdit!")

if __name__ == "__main__":
    main()
