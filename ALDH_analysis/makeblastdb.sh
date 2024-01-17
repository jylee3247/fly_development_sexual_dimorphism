#!/usr/bin/env bash
# Make input FASTA format sequences to blast database
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

INPUT=$1 # Input protein sequence path
OUTPUT=$2 # Output path

# Build BLAST database
makeblastdb -in $INPUT -dbtype prot -parse_seqids -out out