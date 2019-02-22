#! /bin/bash

ASSEMBLER="trinity"

GIT_DIR=$(dirname $(readlink -f ./Snakefile))"/.git"
GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
OUTFILE="transXpress-$ASSEMBLER.stdout.log"
ERRFILE="transXpress-$ASSEMBLER.stderr.log"
echo "$(date)" | tee -a $OUTFILE
echo "transXpress-snakemake now running. git hash: "${GIT_HASH} | tee -a $OUTFILE
echo "Logs are being written to $OUTFILE and $ERRFILE in the current directory" | tee -a $OUTFILE
echo "Try 'lsof $OUTFILE' if you need to get the process id of the snakemake manager" | tee -a $OUTFILE
echo "'tail -f $OUTFILE' will let you see the output of nextflow in real time" | tee -a $OUTFILE 
echo "transXpress-nextflow dropping to background on host "$HOSTNAME"..." | tee -a $OUTFILE

./run_tak.sh 1>$OUTFILE 2>$ERRFILE &
disown
