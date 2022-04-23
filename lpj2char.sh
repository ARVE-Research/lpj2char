#!/bin/bash

datadir=/home/terraces/projects/Emeline_borealforest/lpj_holocene/LPJ-LMfire/output

for infile in `ls $datadir/*.nc`
do
  fname=${infile##*/}

  outfile=${fname%%.*}_zscore.nc
  
  ncgen -k 4 -o $outfile zscore.cdl
  
  echo "working on $outfile"  >&2
  
  ./postprocess $infile $outfile

done
