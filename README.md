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

2. Install [Mambaforge](https://github.com/conda-forge/miniforge#mambaforge)
~~~~
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-Linux-x86_64.sh"
bash Mambaforge-Linux-x86_64.sh
rm Mambaforge-Linux-x86_64.sh
~~~~

3. To ensure correct versions of R packages will be used unset R_LIBS_SITE
~~~~
unset R_LIBS_SITE
~~~~

4. Setup main transXpress conda environment:
~~~~
mamba activate base
mamba create -c conda-forge -c bioconda --name transxpress
mamba activate transxpress
conda config --set channel_priority disabled
mamba env update --file envs/default.yaml
~~~~
* (TIP: if you have a problem updating the default environment try putting python 3.9 into the *defaults.yaml* file)


6. Create a tab-separated file called *samples.txt* in the assembly directory describing where to find your raw read FASTQ files. Create this file with the following contents:
	~~~
	cond_A	cond_A_rep1	A_rep1_left.fq	A_rep1_right.fq
	cond_A	cond_A_rep2	A_rep2_left.fq	A_rep2_right.fq
	cond_B	cond_B_rep1	B_rep1_left.fq	B_rep1_right.fq
	cond_B	cond_B_rep2	B_rep2_left.fq	B_rep2_right.fq
	~~~
    

* You can download reads from SRA with provided script:
    ~~~~
    ./sra_download.sh <SRR0000000> <SRR0000001> <...>
    ~~~~
    where SRR0000000 is an SRA readset ID. E.g., SRR3883773

* See the tests directory for an example of a samples file: [test_samples.txt](./tests/test_samples.txt)

    
7. Take a look at the configuration file *config.yaml* and update as required (you should check: **assembler, targetp, signalp, strand_specific, lineage**).


8. Setup other conda environments (This will take a while):
~~~~
snakemake --use-conda --conda-frontend mamba --conda-create-envs-only --cores 10
~~~~

9. Install SignalP 6.0 (fast):
      * Download SignalP 6.0 fast from https://services.healthtech.dtu.dk/services/SignalP-6.0/ (go to Downloads)
      * Unpack and install signalp:
        ~~~~
         tar zxvf signalp-6.0h.fast.tar.gz
         cd signalp6_fast
         pip install signalp-6-package/
         SIGNALP_DIR=$(python -c "import signalp; import os; print(os.path.dirname(signalp.__file__))" )
         cp -r signalp-6-package/models/* $SIGNALP_DIR/model_weights/
        ~~~~
        (make sure the conda python is used, or use the full path to python from your conda installation)

10. Install TargetP 2.0:
      * Download TargetP 2.0 from https://services.healthtech.dtu.dk/software.php
      * extract the tarball and add path to targetp /bin/ folder to the PATH variable
        ~~~~
         tar zxvf targetp-2.0.Linux.tar.gz
         export PATH=$PATH:`pwd`/targetp-2.0/bin
        ~~~~

## Running transXpress

There are several options on how to run transXpress

Option 1 - use the provided script:
~~~~
./transXpress.sh
~~~~

Option 2 - run snakemake manually with 10 local threads:
~~~~
snakemake --conda-frontend conda --use-conda --cores 10
~~~~

Option 3 - run snakemake manually on an LSF cluster:
~~~~
snakemake --conda-frontend conda --use-conda --latency-wait 60 --jobs 10000 --cluster 'bsub -n {threads} -R "rusage[mem={params.memory}000] span[hosts=1]" -oo {log}.bsub'
~~~~

Option 4 - define a profile and run snakemake with the profile:
~~~~
snakemake --profile profiles/ "$@"
~~~~
You can find example of a simple profile `config.yaml` for slurm in the `profiles` directory or [here](https://github.com/snakemake-profiles/doc)

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

or when using profiles:
~~~
# run only the rules until the multiqc_before_trim rule to check quality of the input data
snakemake multiqc_before_trim --profile profiles/ "$@"
~~~
~~~
# run only the rules to produce samples_trimmed.txt file
snakemake samples_trimmed.txt --profile profiles/ "$@"
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
