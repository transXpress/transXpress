# transXpress-snakemake
transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-nextflow](https://github.com/transXpress/transXpress-nextflow)

## Intro

## Dependencies

Requires
* snakemake 5.4.2+ (install via conda)
* Trinity (install via conda)
* TransDecoder (install via conda)
* BioPython (install via conda)
* samtools (install via conda)
* bowtie2 (install via conda)
* infernal (install via conda)
* HMMER (install via conda)
* kallisto (install via conda)
* [deeploc](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?deeploc)
* [R](https://www.r-project.org)
* [NCBI BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
* basic Linux utitilies: wget, split, awk

## Installation

1. Install [Miniconda3](https://conda.io/en/latest/miniconda.html)
2. Install snakemake and other dependencies:  
```
 conda config --add channels bioconda
 conda config --add channels conda-forge
 conda install snakemake trinity transdecoder biopython samtools bowtie2 infernal hmmer kallisto
```

## Usage

Run locally with 10 parallel threads:
~~~~
snakemake --cores 10 --config samples_file=samples_file.txt
~~~~

Run on an LSF cluster:
~~~~
snakemake --latency-wait 60 --jobs 10000 --cluster 'bsub -n {threads} -R "rusage[mem={params.memory}000] span[hosts=1]" -oo {log}.bsub'
~~~~

