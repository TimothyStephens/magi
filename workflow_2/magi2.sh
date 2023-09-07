#!/usr/bin/env bash

VERSION="0.1"

## Pre-run setup
set -euo pipefail
IFS=$'\n\t'

# Absolute path to this script, e.g. /home/user/bin/foo.sh
SCRIPT=$(readlink -f "$0")
# Absolute path this script is in, thus /home/user/bin
SCRIPTPATH=$(dirname "$SCRIPT")


## Helper functions
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


## Set option envs
NCPUS=8
MIN_DIAMETER=12 # this is the Retro Rules diameter, see https://retrorules.org/doc for details
MAGI_PATH="${SCRIPTPATH}"
FASTA=""
SMILES=""
MZ=""
OUTPUT_DIRECTORY="outout_magi2"


## Useage information
usage() {
echo -e "##
## $(basename ${0}) v${VERSION}
##

Run MAGI2

Usage: 
./$(basename $0) -f pep.fa -m mz.txt
*OR*
./$(basename $0) -f pep.fa -s SMILES.csv

Required:
-f, --fasta     Input fasta file with protein sequences and unique identifiers in the header

-m, --mz        M/Z values (single column, no column header) 
*OR*
-s, --smiles    Compounds in SMILES format (needs a column called 'original_compound')


Optional:
--magi_path     Full path to location where magi is installed (default: ${MAGI_PATH})
-o, --output    Output directory with MAGI2 results (default: ${OUTPUT_DIRECTORY})
-n, --ncpus     Num threads to use for MAGI2 (default: ${NCPUS})
-v, --version   Script version (v${VERSION})
-h, --help      This help message
--debug         Run debug mode
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


## Check required arguments are set
set +eu
if [ -z "${FASTA}" ] || ([ -z "${MZ}" ] && [ -z "${SMILES}" ]);
then
    usage
fi
set -eu



## Check files exist
if [ ! -s "${FASTA}" ]; then log "${FASTA} does not exist!"; exit 1; fi
if [ ! -z "${MZ}" ];
then
if [ ! -s "${MZ}" ]; then log "${MZ} does not exist!"; exit 1; fi
else
if [ ! -s "${SMILES}" ]; then log "${SMILES} does not exist!"; exit 1; fi

fi


## Start running MAGI2

# Create output dir if doesnt exist.
mkdir -p "${OUTPUT_DIRECTORY}" "${OUTPUT_DIRECTORY}/.checkpoint"

# Run "mz_to_InChI.py" if m/z values given not SMILES
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


# Start running MAGI2
CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/compound_to_reaction.done"
if [ ! -e "${CHECKPOINT}" ];
then
  run_cmd "python ${MAGI_PATH}/compound_to_reaction.py --compounds ${COMPOUNDS} --fasta ${FASTA} --diameter ${MIN_DIAMETER} --cpu_count ${NCPUS} --use_precomputed_reactions True --output ${OUTPUT_DIRECTORY}" \
    && touch "${CHECKPOINT}" || exit 1
fi

CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/gene_to_reaction.done"
if [ ! -e "${CHECKPOINT}" ];
then
  run_cmd "python ${MAGI_PATH}/gene_to_reaction.py --not_first_script --output ${OUTPUT_DIRECTORY}" \
    && touch "${CHECKPOINT}" || exit 1
fi

CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/reaction_to_gene.done"
if [ ! -e "${CHECKPOINT}" ];
then
  run_cmd "python ${MAGI_PATH}/reaction_to_gene.py --not_first_script --output ${OUTPUT_DIRECTORY}" \
    && touch "${CHECKPOINT}" || exit 1
fi

CHECKPOINT="${OUTPUT_DIRECTORY}/.checkpoint/scoring.done"
if [ ! -e "${CHECKPOINT}" ];
then
  run_cmd "python ${MAGI_PATH}/scoring.py --not_first_script --output ${OUTPUT_DIRECTORY}" \
    && touch "${CHECKPOINT}" || exit 1
fi


## Done
log "Finished running MAGI2!"


