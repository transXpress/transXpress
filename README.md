![transXpress](logo/Transxpress_Logo_RGB.png)

transXpress: a [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline for streamlined de novo transcriptome assembly and annotation

## Intro

## Dependencies

transXpress requires:
* snakemake 5.4.2+ (install via conda, `envs/default.yaml`)
* fastqc (install via conda, `envs/qc.yaml`)
* multiqc (install via conda, `envs/qc.yaml`)
* trimmomatic (install via conda, `envs/trimmomatic.yaml`)
* Trinity (install via conda, `trinity_utils.yaml`)
* SPAdes (install via conda, `envs/rnaspades.yaml`)
* TransDecoder (install via conda, `transdecoder.yaml`)
* BioPython (install via conda, `envs/default.yaml`)
* samtools (install via conda, `envs/trinity_utils.yaml`)
* bowtie2 (install via conda, `envs/trinity_utils.yaml`)
* infernal (install via conda, `envs/rfam.yaml`)
* HMMER (install via conda, `envs/pfam.yaml`)
* kallisto (install via conda, `envs/trinity_utils.yaml`)
* NCBI BLAST+ (install via conda, `envs/blast.yaml`)
* R (install via conda, `envs/trinity_utils.yaml`)
* edgeR (install via conda, `envs/trinity_utils.yaml`)
* seqkit (install via conda, `envs/default.yaml`, `envs/rnaspades.yaml`)
* wget (install via conda, `envs/default.yaml`)
* sra-tools (install via conda, `envs/default.yaml`)
* tidyverse (required for Trinity, install via conda, `envs/trinity_utils.yaml`)
* python, numpy, pip (install via conda, `envs/default.yaml`)
* busco 4+ (install via conda, `envs/busco.yaml`)
* rsem (install via conda, `envs/trinity_utils.yaml`)
* [SignalP 6.0](https://services.healthtech.dtu.dk/service.php?SignalP)
* [TargetP 2.0](https://services.healthtech.dtu.dk/service.php?TargetP-2.0)
* tmhmm.py (install via pip, `envs/default.yaml`)
* basic Linux utitilies: split, awk, cut, gzip

The conda dependencies are installed in smaller conda environments automatically by transXpress (based on yaml files in the `envs` directory). 

## Installation

1. Checkout the transXpress code into the folder in which you will be performing your assembly:
~~~~
git clone https://github.com/transXpress/transXpress.git
~~~~

2. Install [Miniconda3](https://conda.io/en/latest/miniconda.html)

3. To ensure correct versions of R packages will be used unset R_LIBS_SITE
~~~~
unset R_LIBS_SITE
~~~~

4. Setup main transXpress conda environment:
~~~~
conda create --name transxpress
conda activate transxpress
~~~~

5. Install snakemake and other dependencies in the main transXpress conda environment:  
~~~~
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority false
conda env update --file envs/default.yaml
~~~~

6. Setup other conda environments (This will take a while):
~~~~
snakemake --conda-frontend conda --use-conda --conda-create-envs-only --cores 1
~~~~

7. Install SignalP 6.0 (fast):
      * Download SignalP 6.0 fast from https://services.healthtech.dtu.dk/service.php?DeepLoc-1.0 (go to Downloads)
      * Unpack and install signalp:
        ~~~~
         tar zxvf signalp-6.0g.fast.tar.gz
         cd signalp6_fast
         pip install signalp-6-package/
         SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
         cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
        ~~~~
        (make sure the conda python is used, or use the full path to python from your conda installation)

8. Install TargetP 2.0:
      * Download TargetP 2.0 from https://services.healthtech.dtu.dk/software.php
      * extract the tarball and add path to targetp /bin/ folder to the PATH variable
        ~~~~
         tar zxvf targetp-2.0.Linux.tar.gz
         export PATH=$PATH:`pwd`/targetp-2.0/bin
        ~~~~

## Input

Create a tab-separated file called *samples.txt* with the following contents:
~~~
cond_A    cond_A_rep1    A_rep1_left.fq    A_rep1_right.fq
cond_A    cond_A_rep2    A_rep2_left.fq    A_rep2_right.fq
cond_B    cond_B_rep1    B_rep1_left.fq    B_rep1_right.fq
cond_B    cond_B_rep2    B_rep2_left.fq    B_rep2_right.fq
~~~

Also take a look at the configuration file *config.yaml* and update as required.

You can download reads from SRA with provided script:
~~~~
./sra_download.sh <SRR0000000> <SRR0000001> <...>
~~~~
where SRR0000000 is an SRA readset ID. E.g., SRR3883773
## Running transXpress

Use the provided script:
~~~~
./transXpress.sh
~~~~

Or run snakemake manually with 10 local threads:
~~~~
snakemake --conda-frontend conda --use-conda --cores 10
~~~~

Or run snakemake manually on an LSF cluster:
~~~~
snakemake --conda-frontend conda --use-conda --latency-wait 60 --jobs 10000 --cluster 'bsub -n {threads} -R "rusage[mem={params.memory}000] span[hosts=1]" -oo {log}.bsub'
~~~~

Or define a profile and run snakemake with the profile:
~~~~
snakemake --profile profiles/ "$@"
~~~~
You can find example of a simple profile `config.yaml` for slurm in the `profiles` directory or [here] (https://github.com/snakemake-profiles/doc)

### Running specific steps

You can run specific steps of the pipeline by specifying a rule or a resulting file:
~~~
# run only the rules until the multiqc_before_trim rule to check quality of the input data
./transXpress.sh multiqc_before_trim
~~~

~~~
# run only the rules to produce samples_trimmed.txt file
./transXpress.sh samples_trimmed.txt
~~~

## Running tests
~~~~
cd tests
./run_test.sh
~~~~

## Align reads to the transcriptome assembly and visualize the results in IGV
If you want to align reads to the transcriptome assembly and visualize the results in IGV, you can use the following commands:
~~~~
./transXpress.sh align_reads
./transXpress.sh IGV
~~~~

Then you can load your transcriptome file to IGV: Genomes -> Load Genome from File -> select the file *transcriptome.fasta*

Your sorted .bam files will be in the bowtie_alignments folder: 
~~~~
bowtie_alignments/{sample}.sorted.bam
bowtie_alignments/{sample}.sorted.bam.bai
~~~~
Load them to IGV: File -> Load from File -> select the *bowtie_alignments/{sample}.sorted.bam* files

## The directed acyclic graph (DAG) of the transXpress pipeline execution

![The directed acyclic execution graph](dag.svg )

## Possible problems when executing on cluster systems

### Time limit
Depending on the setup of your cluster you may have to add time option in the transXpress.sh script.
If default time (depends on your cluster setup) for submitted job is not sufficient the job may be cancelled due to time limit.

If this is the case, add time option in submission command in transXpress.sh script.

For example, in case of Slurm change
~~~~
snakemake --conda-frontend conda --use-conda --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "sbatch -o {log}.slurm.out -e {log}.slurm.err -n {threads} --mem {params.memory}GB" "$@"
~~~~
to
~~~~
snakemake --conda-frontend conda --use-conda --latency-wait 60 --restart-times 1 --jobs 10000 --cluster "sbatch -o {log}.slurm.out -e {log}.slurm.err -n {threads} --mem {params.memory} --time=06:00:00" "$@"
~~~~
This sets time limit to 6 hours. You may have to use different time limit based on size of reads used for assembly. 

See https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Computing-Requirements

### Pipeline hangs when cluster cancels the job
It is possible that cluster cancels the job, but pipeline seems to be still running. This can happen because the pipeline does not receive information whether cluster job completed successfully, failed or is still running. You can add `--cluster-status` option and add script which detects the job status. 

See https://snakemake.readthedocs.io/en/stable/tutorial/additional_features.html#using-cluster-status

Alternatively, you can use snakemake [profiles](https://github.com/Snakemake-Profiles/doc) which also contain status checking script. 

See https://snakemake.readthedocs.io/en/v5.1.4/executable.html#profiles 


