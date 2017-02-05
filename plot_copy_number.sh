#!/bin/bash

dir=`pwd`
dir=${dir/FIBROBLASTS/}

Rscript $dir/plot_copy_number.R $1 