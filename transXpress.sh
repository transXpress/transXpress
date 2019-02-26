#! /bin/bash

echo "Running the transXpress pipeline using snakemake"
if [ ! -z `which bsub` ]; then
  snakemake --latency-wait 60 --jobs 10000 --cluster "bsub -oo {log}.bsub -n {threads} -R rusage[mem={params.memory}000] -R span[hosts=1]" --report snakemake_report.html "$@"
else
  snakemake --cores 8 --report snakemake_report.html "$@"
fi


