#!/bin/bash

f=$1
dir=/home/exacloud/lustre1/users/lazar/bin/fastx

name=${f/.fq.gz/}; name=$(basename $name)
gunzip -c $f | \
$dir/fastx_collapser -Q33 > COLLAPSED_SE/$name.uniq.fa
gzip COLLAPSED_SE/$name.uniq.fa

