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
  OUTPUT="${BASENAME}.out"
  OUTPUT_NEW="${OUTPUT}.new"
  COMMAND="${NEWTON_BIN} ${ITER_SCC} -f ${FILE}"
  echo
  echo "Running: ${COMMAND}..."
  time (${COMMAND} >& ${OUTPUT_NEW})
  if [ -f "${OUTPUT}" ]; then
    echo -n "Checking output:"
    diff -q ${OUTPUT} ${OUTPUT_NEW} >& /dev/null
    if [ $? -ne 0 ]; then
      echo " error!"
      echo "Files ${OUTPUT} and ${OUTPUT_NEW} differ!"
    else
      echo " ok."
      rm -f ${OUTPUT_NEW}
    fi
  else
    mv ${OUTPUT_NEW} ${OUTPUT}
  fi
done
