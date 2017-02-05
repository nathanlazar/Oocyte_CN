#!/bin/bash

dir=`pwd`
dir=${dir/FIBROBLASTS/}

Rscript $dir/count_bin_reads.R $1 $2 $3