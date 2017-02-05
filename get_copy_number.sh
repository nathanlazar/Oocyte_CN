#!/bin/bash

dir=`pwd`
dir=${dir/FIBROBLASTS/}

Rscript $dir/get_copy_number.R $1 