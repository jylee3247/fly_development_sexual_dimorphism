#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# Filtering contigs by length
import argparse
import os
from Bio import SeqIO
from glob import glob
from pathlib import Path

# Setup argument parser
parser = argparse.ArgumentParser(description="Filtering out contigs length below 1000.")
parser.add_argument("-i", "--input", type=str, required=True, help="Input contig fasta file.")
parser.add_argument("-o", "--outdir", type=str, required=True, help="Output directory.")
parser.add_argument("-s", "--suffix", type=str, default="_filtered", help="Output prefix")

args = parser.parse_args()

print(f"Input file: {args.input}.\nFiltering out contig length below 1000.")
recs = list(SeqIO.parse(args.input, "fasta")) # Load input fasta file
f_pfx = Path(args.input).stem

# Filter scafolds
new_recs = [rec for rec in recs if len(rec) >= 1000] # Retain contigs with a length of >= 1000

# Save the filtered fasta file
SeqIO.write(new_recs, os.path.join(args.outdir, f_pfx + args.suffix + ".fasta"), "fasta") 

