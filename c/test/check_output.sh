#!/bin/bash

USAGE="Usage: ./check_output.sh <path to newton> \"<options>\" <files>"

NEWTON_BIN="${PWD}/$1"
ITER_SCC=$2
FILES=$3

if [ $# -lt 3 ]; then
  echo "${USAGE}"
  exit 1
fi

for FILE in ${FILES}; do
  BASENAME=$(basename ${FILE})
  SEMIRING=$(echo ${BASENAME} | awk -F _ '{ print $1 }')
  OUTPUT="${BASENAME}.out"
  OUTPUT_NEW="${OUTPUT}.new"
  COMMAND="${NEWTON_BIN} --${SEMIRING} ${ITER_SCC} -f ${FILE}"
  echo
  echo "Running: ${COMMAND}..."
  time (${COMMAND} >& ${OUTPUT_NEW})
done
