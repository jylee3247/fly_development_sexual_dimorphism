#!/usr/bin/env bash
# Multiple sequence alignment using Clustal Omega v1.2.4 
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

INPUT=$1 # Input protein sequence path
OUTPUT=$2 # Output path

# Run Clustal Omega
clustalo-1.2.4 -i $INPUT -o $OUTPUT --distmat-out=${OUTPUT}.msa --full --percent-id