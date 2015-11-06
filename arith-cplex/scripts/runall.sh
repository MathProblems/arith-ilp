#!/bin/bash

## Author: Ashish Sabharwal, AI2, November 2015

## A script to run arithCplex executable on every math problem ILP input file in the
## specified folder. The output for each math problem is saved in a text file with
## the same name as the input file extended with a suffix. The script by default
## expects to find ../arithCplex and ../weights.conf, and launches CPLEX in
## multi-threaded mode.

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
nThreads=8
timelimit=30    # seconds

inputFiles=$datafolder/*$insuffix

time for f in $inputFiles; do
  echo $f
  $arithCplex -t $timelimit --threads $nThreads -s $nSolns --wts $weightsFile $f | /bin/grep -B3 -A1 EXPR > $f$outSuffix
done

