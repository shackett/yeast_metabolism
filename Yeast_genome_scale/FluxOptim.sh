#!/bin/sh

for a_run in 0 1 2 3 4 5 6 7 8 9
do
  qsub -l 1day -cwd -sync n Rscript FluxOptimClus.R runNum=$a_run
done
