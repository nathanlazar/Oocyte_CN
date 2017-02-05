#!/bin/bash

# Usage: make_bins.sh <reads.bam> <reads per bin> <out.bins>

dir=`pwd`
dir=${dir/FIBROBLASTS/}

Rscript $dir/make_bins.R $1 $2 > $3