#!/usr/bin/env bash

# Print all info to log file
exec 1> >(tee "${0}.log.$(date +%s)") 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate ~/miniconda3/envs/magi_2; set -eu
export PYTHONPATH="/home/timothy/miniconda3/envs/magi_2/lib/python3.8/site-packages"

FASTA="test_pep.fa"
MZ="test_mz.txt"

#### Start Script
./../magi2.sh -f "$FASTA" -m "$MZ" --ncpus 48 --nlines 20


