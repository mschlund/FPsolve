#!/bin/bash


USAGE="Usage: ./run_perf.sh <path to newton> <path to test dir> <runs> <options>"

NEWTON_BIN="${PWD}/$1"
GRAMMARS_DIR="${PWD}/$2/grammars"
RUNS=$3
OPTIONS="$4"

if [ ${RUNS} -eq 0 ]; then
  RUNS=1
fi

if [ "x${OPTIONS}" = "x" ]; then
  OPTIONS=""
fi

if [ $# -ne 3 -a $# -ne 4 ]; then
  echo "${USAGE}"
  exit 1
fi

for FILE in "${GRAMMARS_DIR}"/*_perf_*.g; do
  BASENAME=$(basename ${FILE})
  SEMIRING=$(echo ${BASENAME} | awk -F _ '{ print $1 }')
  OUTPUT="${BASENAME}.out"
  COMMAND="${NEWTON_BIN} --${SEMIRING} ${OPTIONS} -f ${FILE}"
  echo
  echo "Running test ${BASENAME}..."
  time (${COMMAND} >& ${OUTPUT})
done
