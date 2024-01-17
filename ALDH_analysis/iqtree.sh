#!/usr/bin/env bash
# Draw phylogenetic tree using IQ-TREE v1.6.12
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

INPUT=$1 # Input MSA path

# Run IQ-TREE
iqtree_path -s $INPUT -nt AUTO

