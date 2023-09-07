#!/usr/bin/env bash

# Print all info to log file
exec 1> "${0}.log.$(date +%s)" 2>&1

#### Pre-run setup
source ~/scripts/script_setup.sh
set +eu; conda activate magi_2; set -eu
export PYTHONPATH="/home/timothy/miniconda3/envs/magi_2/lib/python3.8/site-packages"

FASTA="Porites_Holoproteome_Sequences_for_MAGI.fasta"
MZ="Pos_Metabolite_MAGI.csv"

#### Start Script
/home/timothy/GitHub/magi/workflow_2/magi2.sh -f "$FASTA" -m "$MZ"


