#!/usr/bin/env bash

VERSION="2.0.2"



##
## Pre-run setup
##
set -euo pipefail
IFS=$'\n\t'

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")



##
## Helper functions
##
function run_cmd(){
  local CMD="$@"
  echo -e "[`date`]\tCMD: $CMD"
  eval $CMD
}

function log(){
  echo -e "[`date`]\tLOG: $@"
}

function err(){
  echo -e "[`date`]\tERROR: $@" >&2
}



##
## Set option envs
##
NCPUS=8
NPARTS=8
MIN_DIAMETER=12 # this is the Retro Rules diameter, see https://retrorules.org/doc for details
MAGI_PATH="${SCRIPTPATH}"
FASTA=""
SMILES=""
MZ=""
OUTPUT_DIRECTORY="output_magi2"



##
## Useage information
##
usage() {
echo -e "##
## $(basename ${0}) v${VERSION}
##

Run MAGI2

Usage: 
./$(basename $0) -f pep.fa -p 20 -n 10 -m mz.txt
*OR*
./$(basename $0) -f pep.fa -p 20 -n 10 -s SMILES.csv 

Required:
-f, --fasta     Input fasta file with protein sequences and unique identifiers in the header

-m, --mz        M/Z values (single column, no column header) 
*OR*
-s, --smiles    Compounds in SMILES format (needs a column called 'original_compound')


Optional:
--magi_path     Full path to location where magi is installed (default: ${MAGI_PATH})
-o, --output    Output directory with MAGI2 results (default: ${OUTPUT_DIRECTORY})
-p, --parts     Num parts to split --smiles/--mz file into before running magi (default: ${NPARTS})
-n, --ncpus     Num threads to use for MAGI2 (default: ${NCPUS})
-v, --version   Script version (v${VERSION})
-h, --help      This help message
--debug         Run debug mode

Details:
The --nparts and --ncpus options allows you to control how parallel to run magi 
and how small each parallel piece should be.
NOTE:
  - MAGI can be quite memory hungry at certain stages. Each part (with ~7000 lines)
    can use upto 50GB of memory. So you can split large compound files into lots of
    parts (e.g., '--nparts 24') and run a few of them (e.g., '--ncpus 6') at a time to lower
    overall memory useage. 
" 1>&2
exit 1
}


# See https://stackoverflow.com/questions/192249/how-do-i-parse-command-line-arguments-in-bash

POSITIONAL=()
while [[ $# -gt 0 ]]; do
  key="$1"

  case $key in
    --magi_path)
      MAGI_PATH="$2"
      shift # past argument
      shift # past value
      ;;
    -f|--fasta)
      FASTA="$2"
      shift # past argument
      shift # past value
      ;;
    -m|--mz)
      MZ="$2"
      shift # past argument
      shift # past value
      ;;
    -s|--smiles)
      SMILES="$2"
      shift # past argument
      shift # past value
      ;;
    -o|--output)
      OUTPUT_DIRECTORY="$2"
      shift # past argument
      shift # past value
      ;;
    -p|--nparts)
      NPARTS="$2"
      shift # past argument
      shift # past value
      ;;
    -n|--ncpus)
      NCPUS="$2"
      shift # past argument
      shift # past value
      ;;
    -h|--help)
      usage
      exit 1;
      ;;
    -v|--version)
      echo "v${VERSION}"
      exit 0;
      ;;
    --debug)
      set -x
      shift # past argument
      ;;
    *) # unknown option
      POSITIONAL+=("$1") # save it in an array for later
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters



##
## Starting MAGI2
##
log "MAGI2 (v${VERSION})"



##
## Check required arguments are set
##
set +eu
if [ -z "${FASTA}" ] || ([ -z "${MZ}" ] && [ -z "${SMILES}" ]);
then
  usage
fi
set -eu



##
## Check files exist
##
if [ ! -s "${FASTA}" ]; then log "${FASTA} does not exist!"; exit 1; fi
if [ ! -z "${MZ}" ];
then
  if [ ! -s "${MZ}" ]; then log "${MZ} does not exist!"; exit 1; fi
else
  if [ ! -s "${SMILES}" ]; then log "${SMILES} does not exist!"; exit 1; fi
fi



##
## Create output dir if doesnt exist.
##
mkdir -p "${OUTPUT_DIRECTORY}" "${OUTPUT_DIRECTORY}/.checkpoint"



##
## Run "mz_to_InChI.py" if m/z values given not SMILES
##
if [ ! -z "${MZ}" ];
then
  COMPOUNDS="${OUTPUT_DIRECTORY}/input.smiles.tsv"
  
  CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/mz_to_InChI.done"
  if [ ! -e "${CHECKPOINT}" ];
  then
    run_cmd "python ${MAGI_PATH}/mz_to_SMILES.py -i ${MZ} -o ${COMPOUNDS}" \
      && touch "${CHECKPOINT}" || exit 1
  fi
else
  COMPOUNDS="${SMILES}"
fi



##
## Split input file to run parallel
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/split.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Splitting ${COMPOUNDS} into ${NPARTS} parts"
  run_cmd "${MAGI_PATH}/split.sh -p ${NPARTS} -f ${COMPOUNDS}" \
    && touch "${CHECKPOINT}" || exit 1
  log "  - Done"
fi



##
## MAGI2 - compound_to_reaction.py
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/compound_to_reaction.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Running compound_to_reaction.py on each part"
  while read F;
  do
    N="${F#*.split_*}"
    N="${N%*.*}"
    CHECKPOINT_N="${OUTPUT_DIRECTORY}/.checkpoint/compound_to_reaction.$N.done"
    if [ ! -e "${CHECKPOINT_N}" ];
    then
      echo "python ${MAGI_PATH}/compound_to_reaction.py --compounds ${F} --fasta ${FASTA} --diameter ${MIN_DIAMETER} --cpu_count 1 --use_precomputed_reactions True --output ${F}.magi 1>${F}.compound_to_reaction.log 2>&1 && touch ${CHECKPOINT_N} || exit 1"
    fi
  done < "${COMPOUNDS}.parts.txt" \
    | parallel -j ${NCPUS} -v \
    && touch "${CHECKPOINT}"
  log "  - Done"
fi



##
## MAGI2 - gene_to_reaction.py
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/gene_to_reaction.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Running gene_to_reaction.py on each part"
  while read F;
  do
    N="${F#*.split_*}"
    N="${N%*.*}"
    CHECKPOINT_N="${OUTPUT_DIRECTORY}/.checkpoint/gene_to_reaction.$N.done"
    if [ ! -e "${CHECKPOINT_N}" ];
    then
      echo "python ${MAGI_PATH}/gene_to_reaction.py --not_first_script --output ${F}.magi 1>${F}.gene_to_reaction.log 2>&1 && touch ${CHECKPOINT_N} || exit 1 "
    fi
  done < "${COMPOUNDS}.parts.txt" \
    | parallel -j ${NCPUS} -v \
    && touch "${CHECKPOINT}"
  log "  - Done"
fi



##
## MAGI2 - reaction_to_gene.py
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/reaction_to_gene.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Running reaction_to_gene.py on each part"
  while read F;
  do
    N="${F#*.split_*}"
    N="${N%*.*}"
    CHECKPOINT_N="${OUTPUT_DIRECTORY}/.checkpoint/reaction_to_gene.$N.done"
    if [ ! -e "${CHECKPOINT_N}" ];
    then
      echo "python ${MAGI_PATH}/reaction_to_gene.py --not_first_scrip --output ${F}.magi 1>${F}.reaction_to_gene.log 2>&1 && touch ${CHECKPOINT_N} || exit 1"
    fi
  done < "${COMPOUNDS}.parts.txt" \
    | parallel -j ${NCPUS} -v \
    && touch "${CHECKPOINT}"
  log "  - Done"
fi



##
## MAGI2 - scoring.py
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/scoring.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Running scoring.py on each part"
  while read F;
  do
    N="${F#*.split_*}"
    N="${N%*.*}"
    CHECKPOINT_N="${OUTPUT_DIRECTORY}/.checkpoint/scoring.$N.done"
    if [ ! -e "${CHECKPOINT_N}" ];
    then
      echo "python ${MAGI_PATH}/scoring.py --not_first_script --output ${F}.magi 1>${F}.scoring.log 2>&1 && touch ${CHECKPOINT_N} || exit 1"
    fi
  done < "${COMPOUNDS}.parts.txt" \
    | parallel -j 1 -v \
    && touch "${CHECKPOINT}"
  log "  - Done"
fi



##
## MAGI2 - combine results
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/combine.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Combining magi_compound_results.csv files from each part"
  while read F;
  do
    cat "${F}.magi/magi_compound_results.csv"
  done < "${COMPOUNDS}.parts.txt" \
    | awk 'NR==1 || $1!~"^original_compound"' \
    > "${OUTPUT_DIRECTORY}.magi_compound_results.csv"
  
  log "Combining magi_gene_results.csv files from each part"
  while read F;
  do
    cat "${F}.magi/magi_gene_results.csv"
  done < "${COMPOUNDS}.parts.txt" \
    | awk 'NR==1 || $1!~"^original_compound"' \
    > "${OUTPUT_DIRECTORY}.magi_gene_results.csv"
  
  log "Combining magi_results.csv files from each part"
  while read F;
  do
    cat "${F}.magi/magi_results.csv"
  done < "${COMPOUNDS}.parts.txt" \
    | awk 'NR==1 || $1!~"^original_compound"' \
    > "${OUTPUT_DIRECTORY}.magi_results.csv"
  
  touch "${CHECKPOINT}"
  log "  - Done"
fi



##
## Filter results
##
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/filter.done"
if [ ! -e "${CHECKPOINT}" ];
then
  log "Filtering combined results"
  run_cmd "python ${MAGI_PATH}/filter_results.py -i ${OUTPUT_DIRECTORY}.magi_results.csv -o ${OUTPUT_DIRECTORY}.filtered_magi_results.csv"
  run_cmd "python ${MAGI_PATH}/filter_results.py -i ${OUTPUT_DIRECTORY}.magi_gene_results.csv -o ${OUTPUT_DIRECTORY}.filtered_magi_gene_results.csv"
  run_cmd "python ${MAGI_PATH}/filter_results.py -i ${OUTPUT_DIRECTORY}.magi_compound_results.csv -o ${OUTPUT_DIRECTORY}.filtered_magi_compound_results.csv"
  touch "${CHECKPOINT}"
  log "  - Done"
fi



##
## Done MAGI
##
log "Finished running MAGI2!"



