#!/bin/bash
# Larry Wilhelm.

CA=./runcutadapt.sh
PCA=./post_cutadapt_cleanup.pl

f1=$1
U=$2  # amount to trim from 5' of R1 
D=$3  # amount to trim from 3' of R1
P=$4  # amount to trim from 5' of R2 
Q=$5  # amount to trim from 3' of R2

f2=$f1
f2=${f2/R1/R2};
f1_base=${f1/.fastq.gz/};
f2_base=${f2/.fastq.gz/};
echo "$f1 $f2";
echo "$f1_base $f2_base";
echo "with $U $D $P $Q";

echo "$CA $f1_base R1";
$CA $f1_base R1;

echo "$CA $f2_base R2";
$CA $f2_base R2;

# echo "$PCA $f1_base.trim $f2_base.trim"; 
$PCA $f1_base.trim $f2_base.trim $U $D $P $Q > $f1_base.stats;

rm $f1_base.trim
rm $f2_base.trim

echo "Filtering Stats:" >> $f1_base.report
cat $f1_base.stats >> $f1_base.report
rm $f1_base.stats

mv $f1_base.trim.clean $f1_base.clean.fastq;
gzip $f1_base.clean.fastq;

mv $f2_base.trim.clean $f2_base.clean.fastq;
gzip $f2_base.clean.fastq;

gzip $f1_base.report;
gzip $f2_base.report;

echo "done cleaning $f1"

