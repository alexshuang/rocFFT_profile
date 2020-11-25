#!/bin/sh

CMD=$1
OUT_DIR=${2:-out}
mkdir -p $OUT_DIR

BASIC_PMC="L2CacheHit"
echo "pmc: $BASIC_PMC" > /tmp/input.txt
rocprof -i /tmp/input.txt --hsa-trace --timestamp on -o $OUT_DIR/basic_prof.csv $CMD

cat $OUT_DIR/basic_prof.csv
