#!/bin/bash
# Larry Wilhelm.

CA=./cutadapt-1.2.1/bin/cutadapt
PCA=./post_cutadapt_cleanup_se.pl

f1=$1
U=$2  # amount to trim from 5' of the read 
D=$3  # amount to trim from 3' of the read

f2=$f1
f1_base=${f1/.fq.gz/};

echo "$f1";
echo "$f1_base";
echo "with $U $D";

echo "$CA $f1_base R1";
R1_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#R2_adaptor=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
$CA -a $R1_adaptor -o $f1_base.trim $f1_base.fq.gz > $f1_base.report

echo "$PCA $f1_base.trim"; 
$PCA $f1_base.trim $U $D > $f1_base.stats;

rm $f1_base.trim

echo "Filtering Stats:" >> $f1_base.report
cat $f1_base.stats >> $f1_base.report
rm $f1_base.stats

mv $f1_base.trim.clean $f1_base.clean.fastq;
gzip $f1_base.clean.fastq;

gzip $f1_base.report;

echo "done cleaning $f1"

