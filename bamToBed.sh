#!/bin/bash

/home/users/lazar/bedtools2/bin/bamToBed -i $1 > ${1/.bam/.bed}