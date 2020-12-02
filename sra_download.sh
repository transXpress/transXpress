#!/bin/sh

if [ $# -lt 1 ]; then
  echo Usage:
  echo "$0 <SRR0000000> <SRR0000001> <...>"
  echo "where SRR0000000 is an SRA readset ID. E.g., SRR3883773"
  echo "Single SRA ID or multiple SRA IDs separated by a space is possible."
  exit 1
fi

# Automatically exit if any command returns error
set -e

# This prevents occasional download errors, see https://www.biostars.org/p/161146/ 
vdb-config --restore-defaults

# This sets the temporary folder, in case there is a quota for the default tmp folder
# vdb-config -s /repository/user/main/public/root=$HOME

# 1) The --split-files parameter is necessary, because otherwise F and R reads are
#    combined into single reads of double length
# 2) The --defline-seq '@$sn[_$rn]/$ri' parameter is for processing with Trinity.
#    See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-FAQ#ques_sra_fq_conversion

for DATASET_ID in "$@"; do
  echo Attempting to download: $DATASET_ID
  fastq-dump -$DATASET_ID --defline-seq '@$sn/$ri' --split-files --gzip
done

