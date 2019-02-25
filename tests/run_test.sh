#! /bin/bash

#ASSEMBLER="trinity"

#GIT_DIR=$(dirname $(readlink -f ./Snakefile))"/.git"
#GIT_HASH=$(git --git-dir=${GIT_DIR} log --pretty=format:'%H' -n 1)
#OUTFILE="transXpress-$ASSEMBLER.stdout.log"
#ERRFILE="transXpress-$ASSEMBLER.stderr.log"
#echo "$(date)" | tee -a $OUTFILE
#echo "transXpress-snakemake now running. git hash: "${GIT_HASH} | tee -a $OUTFILE
#echo "Logs are being written to $OUTFILE and $ERRFILE in the current directory" | tee -a $OUTFILE
#echo "Try 'lsof $OUTFILE' if you need to get the process id of the snakemake manager" | tee -a $OUTFILE
#echo "'tail -f $OUTFILE' will let you see the output of nextflow in real time" | tee -a $OUTFILE 
#echo "transXpress-nextflow dropping to background on host "$HOSTNAME"..." | tee -a $OUTFILE

SNAKEFILE="../Snakefile-trinity"

echo "Running the transXpress-trinity pipeline using snakemake"
if [ ! -z `which bsub` ]; then
  snakemake --snakefile $SNAKEFILE --latency-wait 60 --jobs 10000 --cluster "bsub -oo {log}.bsub -n {threads} -R rusage[mem={params.memory}000] -R span[hosts=1]" "$@"
else
  snakemake --snakefile $SNAKEFILE --cores 8 "$@"
fi

echo "Making DAG file describing pipeline execution"
snakemake --snakefile $SNAKEFILE --dag | dot -Tsvg > dag.svg

