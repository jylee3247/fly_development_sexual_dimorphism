#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# SAMtools version 1.16.1
import argparse
import os
import subprocess
import re

def run_samtools(input_path, outdir, n_threads = 4):
    print("Run SAMtools.")
    os.makedirs(outdir, exist_ok=True)


    # Convert SAM to BAM
    out_bam_path = os.path.join(outdir, re.sub("sam$", "bam$", os.path.basename(input_path)))

    cmd = ["samtools", "view", "-b", "--threads", str(n_threads), "-o", out_bam_path, input_path]

    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"SAMtools SAM to BAM convert error: {e.stderr}")

    # SAMtools sorting BAM
    sort_bam_path = os.path.join(outdir, re.sub("bam$", "sorted.bam$", os.path.basename(out_bam_path)))

    cmd = ["samtools", "sort", "--threads", str(n_threads), "-o", sort_bam_path, out_bam_path]
    
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"SAMtools BAM sorting error: {e.stderr}")

# Arguments
parser = argparse.ArgumentParser(description="Run SAMtools.")
parser.add_argument("-i", "--input_path", type=str, required = True, help="Input SAM file path.")
parser.add_argument("-o", "--outdir", type=str, required = True, help="Output directory.")
parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")
args = parser.parse_args()

# Run SAMtools
run_samtools(args.input_path, args.outdir, args.n_threads)