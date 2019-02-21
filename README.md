# transXpress-snakemake
transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-nextflow](https://github.com/transXpress/transXpress-nextflow)

## Intro

## Dependencies

Requires
* snakemake 5.4.2+ (install via conda)
* BioPython (install via conda)
* samtools (install via conda)
* bowtie2 (install via conda)
* infernal (install via conda)
* [R](https://www.r-project.org)
* [NCBI BLAST+](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)
* basic Linux utitilies: wget, split, awk

## Installation

1. Install miniconda: https://conda.io/en/latest/miniconda.html
2. Install snakemake and other dependencies:
```conda install -c bioconda -c conda-forge snakemake biopython samtools bowtie2 infernal```

## Usage

~~~~
snakemake --config samples_file=samples_file.txt
~~~~


