#!/bin/sh

snakemake --latency-wait 60 --jobs 10000 --cluster-sync 'bsub -K -n {threads} -R "rusage[mem={params.memory}000] span[hosts=1]" -oo {log}.bsub' "$@"

