#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# MetaPhlAn4 version 4.0.1
# Blanco-Míguez et al., 2023. doi:10.1038/s41587-023-01688-w
import argparse
import os
import subprocess

def run_metaphlan(
    mpa_inputs,
    out_prefix,
    output_dir,
    analysis_type,
    nproc,
    bowtie2db,
    unclassified_estimation,
    ignore_eukaryotes,
    ignore_archaea,
    input_type
):
    
    print("=====Run MetaPhlAn4=====")
    os.makedirs(output_dir, exist_ok=True)
    
    # Settings
    output_path = os.path.join(output_dir, out_prefix + "_profile.txt")
    bowtie2out = os.path.join(output_dir, out_prefix + ".bowtie2.bz2")
    top_output_dir = os.path.dirname(output_dir)

    cmd = ["metaphlan", mpa_inputs, "--input_type", input_type, "-t", analysis_type, "-o", output_path, "--nproc", str(nproc), "--bowtie2out", bowtie2out, "--bowtie2db", bowtie2db]

    if ignore_eukaryotes:
        cmd += ["--ignore_eukaryotes"]
    if ignore_archaea:
        cmd += ["--ignore_archaea"]
    if unclassified_estimation:
        cmd += ["--unclassified_estimation"]

    print(f"\tCurrent file is {out_prefix}")
    print('\tCommand:', end=' ')
    print(' '.join(cmd), end='\n\n')

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        raise Exception(f"Error occurred while running MetaPhlAn.\n\t\t{e.output}")

# Arguments
parser = argparse.ArgumentParser(description="Run MetaPhlAn4.")
parser.add_argument("inputs", type=str, help="Input file(s).")
parser.add_argument("out_prefix", type=str, help="Output file prefix.")
parser.add_argument("output_dir", type=str, help="Output directory")
parser.add_argument("-at", "--analysis_type", type=str, default="rel_ab", help="Output type.")
parser.add_argument("--nproc", type=int, default=2, help="Number of processing cores.")
parser.add_argument("-db", "--bowtie2db", type=str, default="DB/mpa_vJan21_CHOCOPhlAnSGB", help="Path to the Bowtie2 database.")
parser.add_argument("-it", "--input_type", type=str, default="fastq", help="Type of input files.")
parser.add_argument("--unclassified_estimation", action='store_true', help="Enable unclassified estimation.")
parser.add_argument("--ignore_eukaryotes", action='store_true', help="Ignore eukaryotes in the analysis.")
parser.add_argument( "--ignore_archaea", action='store_true', help="Ignore archaea in the analysis.")

args = parser.parse_args()

# Run MetaPhlAn4
run_metaphlan(args.inputs, args.out_prefix, args.output_dir, args.analysis_type, args.nproc, args.bowtie2db, args.unclassified_estimation, args.ignore_eukaryotes, args.ignore_archaea, args.input_type)