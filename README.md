# transXpress-snakemake
transXpress: a Snakemake pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-nextflow](https://github.com/transXpress/transXpress-nextflow)

## Intro

## Dependencies

Requires
* snakemake 5.4+
* BioPython
* samtools
* R
* infernal
* basic linux utitilies: wget, split

## Installation

1. Install conda:
2. Install snakemake and other dependencies:
  ```conda install -c bioconda -c conda-forge snakemake biopython infernal```

## Usage

~~~~
snakemake --config samples_file=samples_file.txt
~~~~


