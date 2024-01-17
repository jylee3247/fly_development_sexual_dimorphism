#!/usr/bin/env bash
# Regrouping HUMAnN3 output to Pfam IDs.

INPUTDIR=$1 # Directory containing humann output files
OUTDIR=$2 # Output directory

mkdir $OUTDIR

for file in $(ls $INPUTDIR)
do
	echo "Input: ${file}"
	humann_regroup_table --input ${file} --output ${OUTDIR}/${file%.txt}_regrouped.txt --groups uniref90_pfam
done
