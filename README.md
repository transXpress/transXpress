# transXpress-snakemake
transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for rapid de novo transcriptome assembly and annotation

Also see our sister project: [transXpress-nextflow](https://github.com/transXpress/transXpress-nextflow)

## Intro

## Dependencies

transXpress Requires
* snakemake 5.4.2+ (install via conda)
* fastqc (install via conda)
* trimmomatic (install via conda)
* Trinity (install via conda)
* SPAdes (install via conda)
* TransDecoder (install via conda)
* BioPython (install via conda)
* samtools (install via conda)
* bowtie2 (install via conda)
* infernal (install via conda)
* HMMER (install via conda)
* kallisto (install via conda)
* NCBI BLAST+ (install via conda)
* R (install via conda)
* edgeR (install via conda)
* seqkit (install via conda)
* wget (install via conda)
* [deeploc](http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?deeploc)
* basic Linux utitilies: split, awk, cut, gzip

## Installation

1. Install [Miniconda3](https://conda.io/en/latest/miniconda.html)
2. Setup conda environment (optional)
~~~~
conda create --name transxpress
conda activate transxpress
~~~~
3. Install snakemake and other dependencies:  
~~~~
 conda config --add channels bioconda
 conda config --add channels conda-forge
 conda config --add channels r
 conda install snakemake fastqc trimmomatic trinity=2.8.5-0 spades transdecoder biopython samtools bowtie2 infernal hmmer kallisto blast r=3.6.0 bioconductor-edger seqkit wget
~~~~
4. Install deeploc
      * Download deeploc from http://www.cbs.dtu.dk/cgi-bin/nph-sw_request?deeploc
      * Install dependencies: `pip install -r requirements.txt`
      * Install deeploc: `python setup.py install` or locally: `python setup.py install --user`

## Input

Create a tab-separated file called *samples.txt* with the following contents: 
  cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
  cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
  cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
  cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq


## Usage

Run locally with 10 parallel threads:
~~~~
snakemake --cores 10 --config samples_file=samples_file.txt
~~~~

Run on an LSF cluster:
~~~~
snakemake --latency-wait 60 --jobs 10000 --cluster 'bsub -n {threads} -R "rusage[mem={params.memory}000] span[hosts=1]" -oo {log}.bsub'
~~~~

## Running tests
~~~~
cd tests
./run_test.sh
~~~~


## Flow

![The directed acyclic execution DAG of transXpress-snakemake-trinity](./tests/dag.svg )
