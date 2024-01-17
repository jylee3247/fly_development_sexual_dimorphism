#!/usr/bin/env bash
# Run blastp v2.14.1+

QUERY=$1 # Input protein sequence path
DB=$2 # Database path
OUTPUT=$3 # Output path

# Run blastp
blastp -query $QUERY -db $DB -num_threads 4 -max_target_seqs 10 -evalue 1e-2 \
-outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore nident gaps qcovs' -out $OUTPUT