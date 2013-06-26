#!/bin/bash


USAGE="Usage: ./run_perf.sh <path to newton> \"<options>\" <files>"

NEWTON_BIN="${PWD}/$1"
OPTIONS="$2"
FILES="${@:3}"

RUNS=5

if [ $# -lt 3 ]; then
  echo "${USAGE}"
  exit 1
fi

for FILE in ${FILES}; do
  BASENAME=$(basename ${FILE})
  OUTPUT="${BASENAME}.log"
  COMMAND="${NEWTON_BIN} ${OPTIONS} -f ${FILE}"
  if [ -f ${OUTPUT} ]; then
    rm -f ${OUTPUT}
  fi
  echo
  echo "Running ${COMMAND}" | tee -a ${OUTPUT}
  for I in {1..5}; do
    echo "Run ${I}:" | tee -a ${OUTPUT}
    OUT=$(${COMMAND} 2>&1)
    echo "${OUT}" >> ${OUTPUT}
    echo "${OUT}" | grep --color=never 'Solving time'
  done
done
