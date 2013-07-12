#!/bin/sh


echo "$1" #provided number of chunks
for a_chunk in `seq 1 $1`
do
  for a_run in `seq 0 9`
  do
    qsub -l 1day -cwd -sync n Rscript FluxOptimClus.R runNum=$a_run chunk=$a_chunk
  done
done

