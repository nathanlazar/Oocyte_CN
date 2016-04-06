#!/bin/bash

F1=$1
dir=/home/exacloud/lustre1/users/lazar/bin/FastUniq/source

F2=${F1/R1/R2}
name1=${F1/.fq.gz/}; name1=$(basename $name1)
name2=${name1/R1/R2}
gunzip $F1
gunzip $F2

# Running FastUniq
echo ${F1/.gz/}  > COLLAPSED/$name1.txt
echo ${F2/.gz/} >> COLLAPSED/$name1.txt
$dir/fastuniq -i COLLAPSED/$name1.txt -o COLLAPSED/$name1.uniq.fq -p COLLAPSED/$name2.uniq.fq

# Clean up 
rm COLLAPSED/$name1.txt
gzip ${F1/.gz/}
gzip ${F2/.gz/}
gzip COLLAPSED/$name1.uniq.fq
gzip COLLAPSED/$name2.uniq.fq
