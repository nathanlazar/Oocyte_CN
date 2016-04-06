#!/bin/bash

F1=$1
dir=/home/exacloud/lustre1/users/lazar/bin/Trimmomatic-0.35

name=${F1/_R1_001.fq.gz/}; name=${name/.R1.fq.gz}; name=$(basename $name)
F2=${F1/R1/R2}
java -jar $dir/trimmomatic-0.35.jar \
  PE -threads 6 -phred33 $F1 $F2 \
  TRIMMED/$name.R1.trim.fq.gz TRIMMED/$name.R1.unpaired.fq.gz \
  TRIMMED/$name.R2.trim.fq.gz TRIMMED/$name.R2.unpaired.fq.gz \
  ILLUMINACLIP:$dir/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 \
  SLIDINGWINDOW:4:15 MINLEN:36
