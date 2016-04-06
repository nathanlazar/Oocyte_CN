#!/bin/bash

F=$1

/home/exacloud/lustre1/users/lazar/bin/FastQC/fastqc --outdir FASTQC_TRIMMED --threads 6 $1
