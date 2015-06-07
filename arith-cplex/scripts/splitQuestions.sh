#!/bin/bash

if [ $# -ne 2 ]; then
  echo "USAGE: splitQuestions.sh unsplitInputFile outdir"
  exit 1
fi

unsplitInputFile=$1
outdir=$2

prefix="q"
suffix=".txt"

mkdir -p $outdir

#split -a3 -d -l7 $unsplitInputFile $prefix
split -a3 -d -l9 $unsplitInputFile $prefix

for f in $prefix???; do
  #tail -n +2 $f > $outdir/$f.txt
  tail -n +3 $f > $outdir/$f.txt
  /bin/rm $f
done

