#!/bin/bash

F1=$1
dir=/home/exacloud/lustre1/users/lazar/bin/Trimmomatic-0.35

name=${F1/_R1_001.fq.gz/}; name=${name/.R1.fq.gz}; name=$(basename $name)
java -jar $dir/trimmomatic-0.35.jar \
  SE -threads 6 -phred33 $F1 \
  TRIMMED/$name.trim.fq.gz  \
  ILLUMINACLIP:$dir/adapters/TruSeq3-SE.fa:2:30:10 LEADING:5 TRAILING:5 \
  SLIDINGWINDOW:4:15 MINLEN:36

