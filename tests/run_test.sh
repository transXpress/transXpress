#! /bin/bash

SNAKEFILE="../Snakefile"

echo "Running the transXpress-trinity pipeline using snakemake"
snakemake --snakefile $SNAKEFILE --cores 8 "$@"

echo "Making DAG file describing pipeline execution"
snakemake --snakefile $SNAKEFILE --dag | dot -Tsvg > dag.svg

