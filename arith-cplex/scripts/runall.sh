#!/bin/bash

if [ $# -ne 1 ]; then
  echo "USAGE: ./runall.sh datafolder"
  exit 1
fi

datafolder=$1

arithCplex=../arithCplex
weightsFile=../weights.conf

insuffix=".txt"                        
outSuffix=".out"
nSolns=1000
nThreads=4
timelimit=30    # seconds

inputFiles=$datafolder/*$insuffix

time for f in $inputFiles; do
  echo $f
  $arithCplex -t $timelimit --threads $nThreads -s $nSolns --wts $weightsFile $f | /bin/grep -B3 -A1 EXPR > $f$outSuffix
done

