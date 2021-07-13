#! /bin/bash

SNAKEFILE="../Snakefile"

MYCPUS="${PBS_NUM_PPN:-$PBS_NCPUS}"
MYCPUS="${NCPUS:-all}"

MEM="${PBS_RESC_MEM:-$SLURM_MEM_PER_NODE}"
if [ ! -z "$PBS_RESC_MEM" ]; then
    MEM=$ENV{"PBS_RESC_MEM"} / 1024 / 1024 / 1024
elif [ ! -z "$SLURM_MEM_PER_NODE" ]; then
    MEM=$SLURM_MEM_PER_NODE
else
    MEM=`free -g | grep 'Mem:' | awk '{print $2}'`
fi

#  the cmdline args override the settings imposed in $SNAKEFILE
if [ -z "$MYCPUS" ]; then
    echo "Running the transXpress-trinity pipeline using snakemake using all locally available CPUs"
    snakemake --snakefile $SNAKEFILE --cores all --resources mem_gb="$MEM" "$@"
else
    echo "Running the transXpress-trinity pipeline using snakemake using $MYCPUS CPUs"
    snakemake --snakefile $SNAKEFILE --cores $MYCPUS --resources mem_gb="$MEM" "$@"
fi

echo "Making DAG file describing pipeline execution"
snakemake --snakefile $SNAKEFILE --dag | dot -Tsvg > dag.svg

