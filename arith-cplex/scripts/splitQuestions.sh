#!/bin/bash

## Author: Ashish Sabharwal, AI2, November 2015

## A simple script that takes a single file containing multiple math problems in
## the ILP input format and splits it into multiple files in the designated output
## directory, each with one math problem, again in the ILP input format. Output
## files are of the form q012.txt, q349.txt, etc.

if [ $# -ne 2 ]; then
  echo "USAGE: splitQuestions.sh unsplitInputFile outdir"
  exit 1
fi

unsplitInputFile=$1
outdir=$2

# generate input filenames of the form q012.txt, q349.txt, etc.
prefix="q"
suffix=".txt"
numDigitsInFilename=3

# each split has 8 lines of which the top 3 and the bottom 1 are irrelevant
numLinesPerSplit=8
topLinesToSkip=3
botLinesToSkip=1

mkdir -p $outdir

split -d -a $numDigitsInFilename -l $numLinesPerSplit $unsplitInputFile $prefix

for f in $prefix???; do
  tail -n +$topLinesToSkip $f | head -n -$botLinesToSkip > $outdir/$f.txt
  /bin/rm $f
done

