#!/bin/bash

name=$1
read=$2
R1_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
R2_adaptor=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

if [ $read = 'R1' ]
then
	cutadapt -a  $R1_adaptor -o $name.trim $name.fastq.gz > $name.report
fi
	
if [ $read = "R2" ]
then
	cutadapt -a  $R2_adaptor -o $name.trim $name.fastq.gz > $name.report
fi


