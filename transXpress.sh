#! /bin/bash

echo "Running the transXpress pipeline using snakemake"

if [ ! -z `which sbatch` ]; then
  echo "Submitting snakemake jobs to SLURM cluster"
  snakemake --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "sbatch -o {log}.slurm.out -e {log}.slurm.err -n {threads} --mem {params.memory}GB" "$@"
else
if [ ! -z `which bsub` ]; then
  echo "Submitting snakemake jobs to LSF cluster"
  snakemake --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "bsub -oo {log}.bsub -n {threads} -R rusage[mem={params.memory}000] -R span[hosts=1]" "$@"
else
  snakemake --cores 8 "$@"
fi


