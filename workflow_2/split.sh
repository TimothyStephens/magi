#!/usr/bin/env bash

VERSION="0.2"

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
NLINES=5000
FILE=""


## Useage information
usage() {
echo -e "##
## $(basename ${0}) v${VERSION}
##

Split input SMILES into parts based on number of lines in input file

Usage: 
./$(basename $0) -p 24 -f input.smiles.tsv

Required:
-f, --file      Inout SMILES file to split. 
-l, --nlines    No. lines in each output file (default: $NLINES)

Optional:
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
    -f|--file)
      FILE="$2"
      shift # past argument
      shift # past value
      ;;
    -l|--nlines)
      NLINES="$2"
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
if [ -z "${FILE}" ];
then
    usage
fi
set -eu



## Check files exist
if [ ! -s "${FILE}" ]; then log "${FILE} does not exist!"; exit 1; fi


## Start running split
EXT="${FILE##*.}"
PREFIX="${FILE%*.*}.split_"
PART_LIST="${FILE}.parts.txt"

log "split.sh (v${VERSION})"

rm -fr "${PART_LIST}"

# Get number of parts to split into using $NLINES
#  - use round robin approach so each file is roughtly even sized
L=$(tail -n +2 "${FILE}" | wc -l)
NPARTS=$(( $(( L / $NLINES )) + 1 ))
tail -n +2 "${FILE}" | split --numeric-suffixes=1 --suffix-length=4 -n r/$NPARTS - "${PREFIX}"

# Add header to each part
for f in "${PREFIX}"*;
do
    head -n 1 "${FILE}" > tmp_file
    cat "$f" >> tmp_file
    rm "$f"
    mv -f tmp_file "$f.tsv"
    echo "$f.tsv" >> "${PART_LIST}"
done



