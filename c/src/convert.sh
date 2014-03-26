#!/bin/bash

filename=$1
cat $filename | sed -e "s/.*node \[shape = doublecircle\]; //" | sed -e "s/.*init -> //" | sed -e "s/[^0-9]*\([0-9]*\) -> \([0-9]*\) \[label=\"\([^\"]*\)\"\]/\1 \2 \3/" | sed -e "/^ .*/d" | grep -v "digraph" | grep -v "^}" | sed -e "s/\\\n/ /g" | sed -e "s/,/ /g" | sed -e "s/;//g"
