#!/bin/bash

name=$1
R1_adaptor=AGATCGGAAGAGCACACGTCTGAACTCCAGTCA
#R2_adaptor=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT

cutadapt-1.2.1/bin/cutadapt -a  $R1_adaptor -o $name.trim $name.fastq.gz > $name.report





