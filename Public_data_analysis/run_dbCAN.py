#!/usr/bin/env python
# coding: utf-8
##### Gut commensal microbe regulated by host Imd signaling  pathway accelerates sexual maturation in female fruit fly #####

# dbCAN v4.0.0 (https://github.com/linnabrown/run_dbcan)
# Zhang et al., 2018. doi:10.1093/nar/gky418
import argparse
import os
import subprocess

def run_dbcan(input_path, out_dir, out_prefix, db_dir, input_type, cgc_finder, n_threads):

    # Build command
    cmd = ["run_dbcan", 
           input_path, input_type,
           "--out_dir", out_dir, 
           "--out_pre", out_prefix,
           "--db_dir", db_dir, 
           "--dia_cpu", str(n_threads),
           "--hmm_cpu", str(n_threads), 
           "--dbcan_thread", str(n_threads),
           "--tf_cpu", str(n_threads),
           "--stp_cpu", str(n_threads),      
          ]
    
    if cgc_finder:
        cmd.extend(["-c", "cluster"])

    # Print status
    print(f"Start 'run_dbcan':")
    print(f"\tInput fasta: {input_path}")
    print(f"\tOutput files: {out_dir}/{out_prefix}")
    print(f"\tInput type: {input_type}")
    print(f"\tRun command:\n\t\t{' '.join(cmd)}")

    # Run command
    try:
        subprocess.run(cmd, check = True)
    except subprocess.CalledProcessError as e:
        with open(os.path.join(os.path.dirname(out_dir), f"{out_prefix}.error.log"), "w") as f:
            f.writelines(e.stderr.decode())
        raise Exception(f"Error occurred: {e.stderr.decode()}")

# Arguments
parser = argparse.ArgumentParser(description="Run dbCAN4.")
parser.add_argument("input_path", type=str, required=True, help="Path to the input file.")
parser.add_argument("out_dir", type=str, required=True, help="Output directory.")
parser.add_argument("out_prefix", type=str, required=True, help="Output prefix.")
parser.add_argument("--db_dir", type=str, default="dbcan_db/", help="Directory of dbCAN database.")
parser.add_argument("--input_type", type=str, default="meta", help="Type of input.")
parser.add_argument("-c", "--cgc_finder", action='store_true', help="Use CGC Finder.")
parser.add_argument("-t", "--n_threads", type=int, default=4, help="Number of threads.")

args = parser.parse_args()

# run_dbCAN
run_dbcan(args.input_path, args.out_dir, args.out_prefix, args.db_dir, args.input_type, args.cgc_finder, args.n_threads)