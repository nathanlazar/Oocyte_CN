#!/bin/bash

f=$1

out=${f/.bam/.bed}
/home/users/lazar/bedtools2/bin/bamToBed -i $f > $out



