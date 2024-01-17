#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# MEGAHIT version 1.2.9.
# Li et al., 2016. doi:10.1093/bioinformatics/btv033
import argparse
import os
import subprocess

def run_megahit(mode, out_dir, out_prefix, input1, input2, unpaired_input, n_threads):
    print("=====Run Megahit=====")
    os.makedirs(out_dir, exist_ok=True)
    
    if mode == "single":
        print(f"Current mode is 'single'. Ignore 'input1' and 'input2' if provided.")
        input_cmd = ["-r", unpaired_input]
    
    elif mode == "paired":
        print(f"Current mode is 'paired'. Also use 'unpaired_input' if provided.")
        input_cmd = ["-1", input1, "-2", input2]
        input_cmd += ["-r", unpaired_input] if unpaired_input else print("'unpaired_input' not provided.")

    print(f"Outputs will be saved in {out_dir}")
    print("\tStart running...")
    cmd = ["megahit", "-o", out_dir, "--out-prefix", out_prefix]
    cmd += input_cmd

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error occurred while running Megahit.\n\t\t{e.output}")

# Arguments
parser = argparse.ArgumentParser(description="Run Megahit.")
parser.add_argument("mode", type=str, choices=["single", "paired"], help="Choose 'single' or 'paired' according to your data type.")
parser.add_argument("out_dir", type=str, help="Output directory.")
parser.add_argument("out_prefix", type=str, help="Output prefix")
parser.add_argument("-1", "--input1", type=str, default=None, help="Path to the first input file for paired-end mode.")
parser.add_argument("-2", "--input2", type=str, default=None, help="Path to the second input file for paired-end mode.")
parser.add_argument("-u", "--unpaired_input", type=str, default=None, help="Path to the input file for single-end mode. Can be optionally provided for the paired-end mode.")
parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")

args = parser.parse_args()

# Run MEGAHIT
run_megahit(args.mode, args.out_dir, args.out_prefix, args.input1, args.input2, args.unpaired_input, args.n_threads)