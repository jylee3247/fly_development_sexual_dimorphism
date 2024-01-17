#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# Mapping metagenomes to contigs using Bowtie2 version 2.5.1
import argparse
import os
import subprocess

def split_reads(reads):
    """
    Split reads according to their file name.
    """
    
    f_read, r_read, un_read = [], [], []

    for r in reads:
        fn = os.path.basename(r)

        if fn.endswith("_1.fastq.gz"):
            f_read.append(r)

        elif fn.endswith("_2.fastq.gz"):
            r_read.append(r)
            
    un_read = list(set(reads) - set(f_read) - set(r_read))    
    f_read = f_read[0] if f_read else None
    r_read = r_read[0] if r_read else None
    un_read = ",".join(un_read) if un_read else None
    
    return f_read, r_read, un_read

def run_bt2(inputs, bt2db, outpath, n_threads, mode, delete_intermediates):
    print("=====Run Bowtie2=====")
    print(f"\tInput: {inputs}\n\tDB: {bt2db}\n\tMode: {mode}\n\n")

    # Split reads
    f_read, r_read, single_read = split_reads(inputs)

    print(f"Input file(s):")
    print(f"\tForward_read: {f_read}\n\tReverse_read: {r_read}\n\tSingle_read(s): {single_read}\n\tDB: {bt2db}\n")

    if f_read and r_read and not single_read:
        mode = "paired"

    elif not f_read and not r_read and single_read:
        mode = "single"

    elif f_read and r_read and single_read:
        mode = "mixed"

    else:
        raise Exception(f"Failed to define input types ('paired', 'single', 'mixed'). f_fastq: {f_read}, r_fastq: {r_read}, s_fastq(s): {single_read}")

    cmd = ["bowtie2", "-x", bt2db, "-S", outpath, "-p", str(n_threads), mode, "--no-unal"]

    if mode == "paired":
        cmd += ["-1", f_read, "-2", r_read]

    elif mode == "single":
        cmd += ["-U", single_read]

    elif mode == "mixed":
        cmd += ["-1", f_read, "-2", r_read, "-U", single_read]

    # Run bowtie2
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Bowtie2 run error: {e.output}")

# Arguments
parser = argparse.ArgumentParser(description="Run Bowtie2.")
parser.add_argument("-i", "--inputs", nargs='+', required = True, help="Input FASTQ file(s).")
parser.add_argument("-x", "--bt2db", type=str, required = True, help="Path to Bowtie2 database.")
parser.add_argument("-o", "--outpath", type=str, required = True, help="Output file path.")
parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")
parser.add_argument("-m", "--mode", type=str, default="--very-sensitive-local", help="Mode for Bowtie2.")
parser.add_argument("-d", "--delete_intermediates", action='store_true', help="Delete intermediate files.")
args = parser.parse_args()

# Run bowtie2
run_bt2(args.inputs, args.bt2db, args.outpath, args.n_threads, args.mode, args.delete_intermediates)