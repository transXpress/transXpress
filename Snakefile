import os
import sys
import shutil
import re
import csv
import Bio.SeqIO
import subprocess
from snakemake.utils import min_version
from Bio.Seq import Seq

min_version("5.4.1")

configfile: "config.yaml"

for executable in ["samtools", "bowtie2", "kallisto", "deeploc", "targetp", "Trinity", "blastx", "blastp", "makeblastdb", "cmscan", "hmmscan", "fastqc", "rnaspades.py", "seqkit", "R"]:
    if not shutil.which(executable):
        sys.stderr.write("Warning: Cannot find %s in your PATH%s" % (executable, os.path.sep))

TRINITY_EXECUTABLE_PATH=shutil.which("Trinity")
TRINITY_HOME=os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))) ##Have to resolve the symbolic link that conda makes

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis

# These rules don't need to be sent to a cluster
localrules: all, clean, trimmomatic_split, trimmomatic_merge, samples_yaml_conversion, trinity_butterfly_split, transcriptome_copy

rule all:
  """
  List of target files of the transxpress pipeline.
  """
  input:
    "samples_trimmed.txt",
    "transcriptome.fasta",
    "transcriptome.pep",
    "transcriptome_stats.txt",
    "ExN50_plot.pdf",
    "transcriptome_annotated.fasta",
    "transcriptome_annotated.pep",
    "transcriptome_TPM_blast.csv",
    "fastqc",
    "multiqc",
    "edgeR_trans",
    "busco",
    "busco_report.txt"


rule clean:
  """
  Removes files from the transxpress directory when specifically called by the user.
  """

  shell:
    """
    if [ -f samples_trimmed.txt ]; then
      cut -f 2 < samples_trimmed.txt | xargs --no-run-if-empty rm -rf
    fi
    rm -rf trinity_* rnaspades_* tmp* log TMHMM* kallisto* transcriptome* pipeliner* annotation* transdecoder* trimmomatic* samples_trimmed* ExN50_plot.pdf multiqc fastqc edgeR_trans
    """


rule fastqc:
  """
  Runs fastQC on individual input reads files.
  """
  input:
    samples=config["samples_file"]
  output:
    directory("fastqc")
  log:
    "logs/fastqc.log"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # run fastqc on input files
    # making directory for the fastqc summary files
    mkdir {output} &> {log}
    FILES=$(awk '{{ printf "%s\\n%s\\n", $3,$4}}' {input}) &>> {log}
    # running fasqc on the files
    for file in $FILES; do fastqc -f fastq -o {output} $file; done &>> {log}
    """


rule multiqc:
  """
  Creates multiQC report from individual fastQC reports.
  """
  input:
    "fastqc"
  output:
    directory("multiqc")
  log:
    "logs/multiqc.log"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    mkdir {output} &> {log}
    multiqc -o {output} {input} &>> {log}
    """


checkpoint trimmomatic_split:
  """
  Splits file with information about all reads files into separate files
  so they can be processed with trimmomatic in parallel.
  """   
  input:
    samples=config["samples_file"]
  output:
    directory("trimmomatic")
  log:
    "logs/trimmomatic_split.log"
  shell:
    """
    mkdir -p {output} &> {log}
    split -d -l 1 {input} {output}/sample_ &>> {log}
    """


rule trimmomatic_parallel:
  """
  Processes individual read files with trimmomatic.
  """
  input:
    "trimmomatic/sample_{job_index}"
  output:
    "trimmomatic/completed_{job_index}"
  log:
    "logs/trimmomatic_parallel{job_index}.log"
  params:
    memory="10"
  threads:
    4
  shell:
    """
    # note: trimmomatic can use gzipped files directly
    read SAMPLE REPLICATE F_READS R_READS < {input}
    # If the sample line is empty, ignore it
    if [ -z "$REPLICATE" ]; then
      touch {output}
      exit 0
    fi
    if [ ! -z "$R_READS" ]; then
      trimmomatic PE -threads {threads} $F_READS $R_READS trimmomatic/{wildcards[job_index]}.R1-P.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R1-U.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R2-P.qtrim.fastq.gz trimmomatic/{wildcards[job_index]}.R2-U.qtrim.fastq.gz {config[trimmomatic_parameters]} &> {log}
      echo $SAMPLE	$REPLICATE	trimmomatic/{wildcards[job_index]}.R1-P.qtrim.fastq.gz	trimmomatic/{wildcards[job_index]}.R2-P.qtrim.fastq.gz > {output} 2>> {log}
    else
      trimmomatic SE -threads {threads} $F_READS trimmomatic/{wildcards[job_index]}.U.qtrim.fastq.gz {config[trimmomatic_parameters]} &> {log}
      echo $SAMPLE      $REPLICATE      trimmomatic/{wildcards[job_index]}.U.qtrim.fastq.gz > {output} 2>> {log}
    fi
    """


def trimmomatic_completed_parallel_jobs(wildcards):
  """
  Returns names of files with information about files processed with trimmomatic.
  """
  parallel_dir = checkpoints.trimmomatic_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "sample_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids


rule trimmomatic_merge:
  """
  Creates a file with information about all reads files processed with trimmomatic.
  """
  input:
    trimmomatic_completed_parallel_jobs
  output:
    samples_trimmed="samples_trimmed.txt"
  log:
    "logs/trimmomatic_merge.log"
  shell:
    """
    cat {input} > {output.samples_trimmed} 2> {log}
    """


rule samples_yaml_conversion:
  """
  Converts the file with information about all trimmed reads files into a yaml
  format used by rnaSPAdes.
  """
  input:
    samples_txt="samples_trimmed.txt"
  output:
    samples_yaml="samples_trimmed.yaml"
  log:
    "logs/samples_yaml_conversion.log"
  run:
    import os
    import os.path
    import pprint

    sample_list = []
    with open(input["samples_txt"], "r") as input_handle:
      for line in input_handle:
        row = re.split("[\t ]+", line)
        if (len(row) > 3): # paired reads
          paired_dict = {}
          paired_dict['orientation'] = 'fr'
          paired_dict['type'] = 'paired-end'
          f_reads = row[2].strip()
          r_reads = row[3].strip()
          assert os.path.isfile(f_reads)
          assert os.path.isfile(r_reads)
          paired_dict['left reads'] = [f_reads]
          paired_dict['right reads'] = [r_reads]
          ##Depends on where the unpaired reads are written to
          ##Presumably could get from sample txt as well
          ##Note: Commented out, to make it consistent with Trinity. May change in future
          ##unpaired_reads_f = f[:-16]+"U.qtrim.fastq.gz"
          ##unpaired_reads_r = r[:-16]+"U.qtrim.fastq.gz"
          ##assert os.path.isfile(unpaired_reads_f)
          ##assert os.path.isfile(unpaired_reads_r)
          ##unpaired_reads = [unpaired_reads_f] + [unpaired_reads_r]
          ##if len(unpaired_reads) > 0:
          ##    sample_dict['single reads'] = unpaired_reads
          sample_list.append(paired_dict)

        if (len(row) == 3): # unpaired reads
          unpaired_dict = {}
          unpaired_dict['type'] = 'single'
          u_reads = row[2].strip()
          assert os.path.isfile(u_reads)
          unpaired_dict['single reads'] = [u_reads]
          sample_list.append(unpaired_dict)

    with open(output["samples_yaml"], "w") as output_handle:
      output_handle.write(pprint.pformat(sample_list))


rule trinity_inchworm_chrysalis:
  """
  Runs first two stages of Trinity assembly (Inchworm and Chrysalis).
  """
  input:
    samples="samples_trimmed.txt",
  output:
    "trinity_out_dir/recursive_trinity.cmds"
  log:
    "logs/trinity_inchworm_chrysalis.log"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --samples_file {input} {config[trinity_parameters]} {config[strand_specific]} &> {log}
    """


checkpoint trinity_butterfly_split:
  """
  Preparation for last stage of Trinity (Butterfly) parallelization by splitting
  the commands into independent parts.
  """
  input:
    "trinity_out_dir/recursive_trinity.cmds"
  output:
    directory("trinity_out_dir/parallel_jobs")
  log:
    "logs/trinity_split.log"
  shell:
    """
    mkdir -p {output} &> {log}
    # Note: it is important to use -d and not --numeric-suffixes, see https://github.com/transXpress/transXpress-snakemake/issues/12
    # Note #2: we split into 1000 chunks to avoid running over the command line limit when too many parallel jobs are created, see https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin 
    split -n l/1000 -e -d {input} {output}/job_ &>> {log}
    """


rule trinity_butterfly_parallel:
  """
  Runs Trinity Butterfly commands (which were split into parts) in parallel.
  """
  input:
    "trinity_out_dir/parallel_jobs/job_{job_index}"
  output:
    "trinity_out_dir/parallel_jobs/completed_{job_index}"
  log:
    "logs/trinity_parallel{job_index}.log"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    bash {input} &> {log}
    cp -p {input} {output} &>> {log}
    """


def trinity_completed_parallel_jobs(wildcards):
  """
  Returns filenames of the files processed in parallel.
  """
  parallel_dir = checkpoints.trinity_butterfly_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "job_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trinity_butterfly_parallel_merge:
  input:
    jobs=trinity_completed_parallel_jobs,
    cmds="trinity_out_dir/recursive_trinity.cmds"
  output:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed"
  log:
    "logs/trinity_butterfly_parallel_merge.log"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    # Can crash if there are too many parallel jobs
    # See https://bitbucket.org/snakemake/snakemake/issues/878/errno-7-argument-list-too-long-path-to-bin 
    cat {input.jobs} > {output.cmds_completed} 2> {log}
    """


rule trinity_final:
  """
  Runs the final Trinity assembly.
  """
  input:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed",
    samples="samples_trimmed.txt"
  output:
    transcriptome="trinity_out_dir/Trinity.fasta",
    gene_trans_map="trinity_out_dir/Trinity.fasta.gene_trans_map"
  log:
    "logs/trinity_final.log"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    Trinity --max_memory {params.memory}G --CPU {threads} --samples_file {input.samples} {config[trinity_parameters]} {config[strand_specific]} &>> {log}
    """


rule rnaspades:
  """
  Runs rnaSPAdes assembly.
  """
  input:
    samples="samples_trimmed.yaml",
  output:
    transcriptome="rnaspades_out/transcripts.fasta",
    ##config["annotated_fasta_prefix"]+"_spades_out/K"+THEKMER+"/assembly_graph_with_scaffolds.gfa"
    gene_trans_map="rnaspades_out/transcripts.gene_trans_map"
  log:
    "logs/rnaspades.log"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    ##TODO = kmer shouldn't be fixed, & should be configurable from the beginning of the script (did run into some bugs with the auto parameter)
    rnaspades.py --dataset {input.samples} -t {threads} -m {params.memory} -o rnaspades_out --only-assembler -k 47 &> {log}
    ##Make a fake gene_trans_map file
    seqkit seq -n {output.transcriptome} | while read id ; do echo -e "$id\\t$id" ; done > {output.gene_trans_map} 2>> {log}
    """
 

rule transcriptome_copy:
  """
  Copies the assembled transcriptome to the transxpress directory.
  """
  input:
    transcriptome=rules.rnaspades.output.transcriptome if config["assembler"]=="rnaspades" else rules.trinity_final.output.transcriptome,
    gene_trans_map=rules.rnaspades.output.gene_trans_map if config["assembler"]=="rnaspades" else rules.trinity_final.output.gene_trans_map,
  output:
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  shell:
    """
    cp -p {input.transcriptome} {output.transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.gene_trans_map} &>> {log}
    """


rule trinity_stats:
  """
  Runs Trinity script to get statistics about the assembled transcriptome
  (number of transcripts, genes, GC content, EXN50).
  """
  input:
    transcriptome="transcriptome.fasta",
    expression="transcriptome_expression_isoform.tsv"
  output:
    stats="transcriptome_stats.txt",
    exN50="transcriptome_exN50.tsv",
    exN50plot="ExN50_plot.pdf"
  log:
    "logs/trinity_exN50.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    {TRINITY_HOME}/util/TrinityStats.pl {input.transcriptome} > {output.stats} 2> {log}
    {TRINITY_HOME}/util/misc/contig_ExN50_statistic.pl {input.expression} {input.transcriptome} > {output.exN50} 2>> {log}
    {TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript {output.exN50} &>> {log}
    """

rule busco:
  """
  Runs BUSCO to assess the completeness of the transcriptome.
  """
  input:
    transcriptome="transcriptome.fasta"
  output:
    out_directory=directory("busco"),
    report="busco_report.txt"
  log:
    "logs/busco.log"
  params:
    memory="2"
  threads:
    4
  shell:
    """
    lineage={config[lineage]} &> {log}
    if [ -z "$lineage"] &>> {log}
    then &>> {log}
      busco -m transcriptome -i {input.transcriptome} -o {output.out_directory} --auto-lineage -c {threads} &>> {log}
    else &>> {log}
      busco -m transcriptome -i {input.transcriptome} -o {output.out_directory} -l $lineage -c {threads} &>> {log}
    fi &>> {log}

    cp $(ls busco/short_summary*.txt) {output.report} &>> {log}
    """

rule transdecoder_longorfs:
  """
  Runs first stage of Transdecoder extracting the long open reading frames.
  """
  input:
    transcriptome="transcriptome.fasta",
  output:
    orfs="transcriptome.orfs"
  log:
    "logs/transdecoder_longorfs.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    rm -rf {input.transcriptome}.transdecoder_dir &> {log}
    TransDecoder.LongOrfs -t {input.transcriptome} --output_dir transdecoder &>> {log} 
    cp -p transdecoder/longest_orfs.pep {output.orfs} &>> {log}
    """


rule transdecoder_predict:
  """
  Runs second stage of Transdecoder predicting the likely coding regions based
  on Pfam and SwissProt hits.
  """
  input:
    transcriptome="transcriptome.fasta",
    pfam="annotations/pfam_orfs.out",
    blastp="annotations/sprotblastp_orfs.out"
  output:
    "transcriptome.pep"
  log:
    "logs/transdecoder_predict.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TransDecoder.Predict -t {input.transcriptome} --output_dir transdecoder --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} &> {log}
    cp -p {input.transcriptome}.transdecoder.pep {output} &>> {log}
    """


rule align_reads:
  """
  Not run by default.
  Runs Trinity script aligning the reads to the transcriptome using Bowtie2
  and estimating the abundance with RSEM.
  """
  input:
    samples=config["samples_file"],
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  output:
    "RSEM_out.gene.TMM.EXPR.matrix"
  log:
    "logs/bowtie2.log"
  params:
    memory="100"
  threads:
    16
  shell:
    """
    {TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts {input[1]} {config[strand_specific]} --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map {input[2]} &> {log}
    # samtools sort -l 9 -@ {threads} -T  ??? ???/bowtie2.bam > bowtie2.sorted.bam
    # samtools index bowtie2.sorted.bam
    """


rule trinity_DE:
  """
  Runs Trinity script to perform differential expression analysis using edgeR.
  """
  input:
    samples=config["samples_file"],
    expression="kallisto.gene.counts.matrix"
  output:
    directory("edgeR_trans")
  log:
    "logs/trinity_DE.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    num_replicates=`awk '{{print $2}}' {input.samples} | sort | uniq | wc -l` &> {log}
    num_samples=`awk '{{print $1}}' {input.samples} | sort | uniq | wc -l` &>> {log}
    num_replicates_minus_samples=$((num_replicates - num_samples)) &>> {log}
    if [ $num_replicates_minus_samples -gt 1 ] &>> {log}
    then &>> {log}
	    {TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} &>> {log}
    else &>> {log}
	    echo "No biological replicates to run proper differential expression analysis, last-resorting to edgeR with --dispersion 0.1" &>> {log}
      {TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} --dispersion {config[dispersion]} &>> {log}
    fi &>> {log}
    """


# a more elegant way to do this is:
# seqkit seq -i {input} | seqkit replace -s -p "\*" -r "" | seqkit split -f -s 1 -O {output}
# however, there seems to be a problem in seqkit: https://github.com/shenwei356/seqkit/issues/65
checkpoint fasta_split_fasta:
  """
  Splits the transcriptome.fasta into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "transcriptome.fasta"
  output:
    directory("annotations/chunks_fasta")
  log:
    "logs/fasta_split.log"
  params:
    memory="2"
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % 1000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 1000) + 1);
          filename = os.path.join(output[0], fileindex + ".fasta")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

checkpoint fasta_split_orfs:
  """
  Splits the transcriptome.orfs into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "transcriptome.orfs"
  output:
    directory("annotations/chunks_orfs")
  log:
    "logs/orfs_split.log"
  params:
    memory="2"
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % 1000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 1000) + 1);
          filename = os.path.join(output[0], fileindex + ".orfs")
          output_handle = open(filename, "w");
        
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        
checkpoint fasta_split_pep:
  """
  Splits the transcriptome.pep into smaller files so
  they can be processed (annotated) in parallel.
  """
  input:
    "transcriptome.pep"
  output:
    directory("annotations/chunks_pep")
  log:
    "logs/pep_split.log"
  params:
    memory="2"
  threads:
    1
  run:
    os.makedirs(output[0], exist_ok=True)
    with open(input[0], "r") as input_handle:
      output_handle = None
      count = 0
      for record in Bio.SeqIO.parse(input_handle, "fasta"):
        if (count % 1000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 1000) + 1);
          filename = os.path.join(output[0], fileindex + ".pep")
          output_handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        #if wildcards["extension"] == "pep":
        record.seq = record.seq.strip("*")
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

def parallel_annotation_tasks_fasta(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_fasta.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.fasta" )).index
  completed_files = expand("annotations/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_tasks_orfs(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_orfs.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.orfs")).index
  completed_files = expand("annotations/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_tasks_pep(wildcards):
  """
  Returns filenames of files which were annotated in parallel.
  """
  parallel_dir = checkpoints.fasta_split_pep.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.pep")).index
  completed_files = expand("annotations/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


rule annotation_merge_fasta:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_fasta
  output:
    "annotations/{task}_fasta.out"
  log:
    "logs/{task}_fasta_merge.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """

rule annotation_merge_orfs:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_orfs
  output:
    "annotations/{task}_orfs.out"
  log:
    "logs/{task}_orfs_merge.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule annotation_merge_pep:
  """
  Merges files that were annotated in parallel into a single file.
  """
  input:
    parallel_annotation_tasks_pep
  output:
    "annotations/{task}_pep.out"
  log:
    "logs/{task}_pep_merge.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule rfam_parallel:
  """
  Runs cmscan on Rfam database on smaller nucleotide files in parallel.
  """
  input:
    fasta="annotations/chunks_fasta/{index}.fasta",
    db="db/Rfam.cm"
  output:
    "annotations/rfam/{index}.out"
  log:
    "logs/rfam_{index}.log"
  params:
    memory="8" # increased memory from 2 to 8 becuase it was not sufficient
  threads:
    2
  shell:
    """
    cmscan -E {config[e_value_threshold]} --rfam --cpu {threads} --tblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule pfam_parallel:
  """
  Runs hmmscan on Pfam database on smaller protein files in parallel.
  """
  input:
    fasta="annotations/chunks_orfs/{index}.orfs",
    db="db/Pfam-A.hmm"
  output:
    "annotations/pfam/{index}.out"
  log:
    "logs/pfam_{index}.log"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    # Transdecoder requires --domtblout output
    hmmscan -E {config[e_value_threshold]} --cpu {threads} --domtblout {output} {input[db]} {input[fasta]} &> {log}
    """


rule sprot_blastp_parallel:
  """
  Runs blast search on SwissProt database on smaller protein files in parallel.
  """
  input:
    fasta="annotations/chunks_orfs/{index}.orfs",
    db="db/uniprot_sprot.fasta"
  output:
    "annotations/sprotblastp/{index}.out"
  log:
    "logs/sprotblastp{index}.log"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastp -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule sprot_blastx_parallel:
  """
  Runs blast search on SwissProt database on smaller nucleotide files in parallel.
  """
  input:
    fasta="annotations/chunks_fasta/{index}.fasta",
    db="db/uniprot_sprot.fasta"
  output:
    "annotations/sprotblastx/{index}.out"
  log:
    "logs/sprotblastx_{index}.log"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastx -query {input[fasta]} -db {input[db]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule tmhmm_parallel:
  """
  Runs tmhmm.py on smaller protein files in parallel to predict transmembrane
  topology.
  """
  input:
    "annotations/chunks_pep/{index}.pep"
  output:
    "annotations/tmhmm/{index}.out"
  log:
    "logs/tmhmm_{index}.log"
  params:
    memory="2"
  threads:
    1
  run:
    output_dir = os.path.join("annotations", "tmhmm")

    # create output directory
    os.makedirs(output_dir, exist_ok=True)

    # open the input file, output file and log file
    with open(input[0], "r") as input_handle, open(output[0], "a+") as output_handle, open(log[0], "a+") as log_handle:

      # iterate through individual sequences in input file, tmhmm.py can be executed with just one sequence at a time
      for record in Bio.SeqIO.parse(input_handle, "fasta"):

        # individual fasta headers must have an ID and description, separated by a space, otherwise error occurs
        # see: https://github.com/dansondergaard/tmhmm.py/issues/16
        # adding artificial "description" in header
        record.description = " x"

        # tmhmm.py does not handle 'X' in the sequence (happens when using RNASpades as an assembler and 'N' occurs in the transcript sequence)
        # see: https://github.com/dansondergaard/tmhmm.py/issues/9
        # removing 'X' from the sequence
        record.seq = Seq(re.sub('X',"",str(record.seq)))

        # saving this fasta file
        fasta_file =  os.path.join(output_dir, record.id + ".fasta")
        Bio.SeqIO.write(record, fasta_file, "fasta")

        # executing the tmhmm.py on created fasta
        cmd_tmhmm = ['tmhmm', '-f', fasta_file]
        subprocess.run(cmd_tmhmm, stdout = log_handle, stderr = log_handle)

        # 3 files are created afterwards
        # record.id.summary     record.id.annotation    record.id.plot

        summary_file = record.id + ".summary"
        annotation_file = record.id + ".annotation"
        plot_file = record.id + ".plot"

        # creating shorter output from summary file
        with open(summary_file, "r") as inp:
          final_topology = ''
          num_of_helices = 0

          for line in inp:
            # split the line of format "start end topology" into list [start, end, topology]
            split_line = re.split(r'\s+', line, 2)
            start = split_line[0].strip()
            end = split_line[1].strip()
            topology = split_line[2].strip()

            # update final string based on the topology
            if (topology == 'inside'):
              final_topology += 'i'
            if (topology == 'outside'):
              final_topology += 'o'
            if ('transmembrane helix' in topology):
              final_topology += start + '-' + end
              num_of_helices += 1

          output_handle.write(record.id + "\tPredHel=" + str(num_of_helices) + "\tTopology=" + final_topology + "\n")
        os.remove(summary_file)
        os.remove(annotation_file)
        os.remove(plot_file)


rule deeploc_parallel:
  """
  Runs Deeploc on smaller protein files in parallel to predict their subcellular
  localization.
  """
  input:
    "annotations/chunks_pep/{index}.pep"
  output:
    "annotations/deeploc/{index}.out"
  log:
    "logs/deeploc_{index}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    # Required for deeploc in some installations
    # See https://github.com/Theano/Theano/issues/6568
    export MKL_THREADING_LAYER=GNU
    deeploc -f {input} -o {output} &> {log}
    mv {output}.txt {output} &>> {log}
    """

rule targetp_parallel:
  """
  Runs TargetP on smaller protein files to predict presence and type of
  targeting peptide.
  """
  input:
    "annotations/chunks_pep/{index}.pep"
  output:
    "annotations/targetp/{index}.out"
  log:
    "logs/targetp_{index}.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    targetp -fasta {input} -format short -org {config[targetp]} -prefix {wildcards.index} &> {log}
    mv {wildcards.index}_summary.targetp2 {output} &>> {log}
    """


rule kallisto:
  """
  Runs Trinity script to perform transcript expression quantification 
  using Kallisto.
  """
  input:
    samples="samples_trimmed.txt",
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  output:
    "transcriptome_expression_isoform.tsv",
    "transcriptome_expression_gene.tsv",
    "kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  params:
    memory="8" # increased memory from 2 to 8 since it was not sufficient
  threads:
    8
  shell:
    """
    {TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} {config[strand_specific]} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
    {TRINITY_HOME}/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input.gene_trans_map} */abundance.tsv &>> {log}
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp -p kallisto.isoform.TMM.EXPR.matrix {output[0]} &>> {log}
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp -p kallisto.isoform.TPM.not_cross_norm {output[0]} &>> {log}
    else
      echo Neither kallisto.isoform.TMM.EXPR.matrix or kallisto.isoform.TPM.not_cross_norm were produced
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp -p kallisto.gene.TMM.EXPR.matrix {output[1]} &>> {log}
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp -p kallisto.gene.TPM.not_cross_norm {output[1]} &>> {log}
    else
      echo Neither kallisto.gene.TMM.EXPR.matrix or kallisto.gene.TPM.not_cross_norm were produced
      exit 1
    fi
    """

 
rule download_sprot:
  """
  Downloads and prepares SwissProt database.
  """
  output:
    "db/uniprot_sprot.fasta"
  log:
    "logs/download_sprot.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" &> {log}
    gunzip db/uniprot_sprot.fasta.gz &>> {log}
    makeblastdb -in db/uniprot_sprot.fasta -dbtype prot &>> {log}
    """


rule download_pfam:
  """
  Downloads and prepares Pfam database.
  """
  output:
    "db/Pfam-A.hmm"
  log:
    "logs/download_pfam.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" &> {log}
    gunzip db/Pfam-A.hmm.gz &>> {log}
    hmmpress db/Pfam-A.hmm &>> {log}
    """


rule download_rfam:
  """
  Downloads and prepares Rfam database.
  """
  output:
    "db/Rfam.cm"
  log:
    "logs/download_rfam.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/Rfam/CURRENT/Rfam.cm.gz" &> {log}
    gunzip db/Rfam.cm.gz &>> {log}
    cmpress db/Rfam.cm &>> {log}    
    """


rule download_eggnog:
  """
  For future development.
  Downloads eggnog database.
  """
  output:
    "db/NOG.annotations.tsv"
  log:
    "logs/download_eggnog.log"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "http://eggnogdb.embl.de/download/latest/data/NOG/NOG.annotations.tsv.gz" &> {log}    
    gunzip db/NOG.annotations.tsv.gz &>> {log}
    """


rule annotated_fasta:
  """
  Puts annotations in the headers of transcripts/proteins in the 
  .fasta/.pep transcriptome files.
  """
  input:
    transcriptome="transcriptome.fasta",
    proteome="transcriptome.pep",
    expression="transcriptome_expression_isoform.tsv",
    blastx_results="annotations/sprotblastx_fasta.out",
    rfam_results="annotations/rfam_fasta.out",
    blastp_results="annotations/sprotblastp_orfs.out",
    pfam_results="annotations/pfam_orfs.out",
    tmhmm_results="annotations/tmhmm_pep.out",
    deeploc_results="annotations/deeploc_pep.out",
    targetp_results="annotations/targetp_pep.out"
  output:
    transcriptome_annotated="transcriptome_annotated.fasta",
    proteome_annotated="transcriptome_annotated.pep",
    tpm_blast_table="transcriptome_TPM_blast.csv"
  log:
    "logs/annotated_fasta.log"
  params:
    memory="4"
  threads:
    1
  run:
    # Open log file
    with open(log[0], "w") as log_handle:
      ## Annotation map: transcript id -> description
      expression_annotations = {}
      blastx_annotations = {}
      blastp_annotations = {}
      pfam_annotations = {}
      rfam_annotations = {}
      tmhmm_annotations = {}
      deeploc_annotations = {}
      targetp_annotations = {}
  
      ## Load kallisto results
      print ("Loading expression values from", input["expression"], file=log_handle)
      with open(input["expression"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter='\t')
        columns = next(csv_reader)
        for row in csv_reader:
          expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
          for i in range(2, len(columns)):
            expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

      ## Load blastx results
      print ("Loading blastx results from", input["blastx_results"], file=log_handle)
      with open(input["blastx_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastx_annotations[row[0]] = row[12] + " E=" + str(row[10])
  
      ## Load blastp results
      print ("Loading blastp results from", input["blastp_results"], file=log_handle)
      with open(input["blastp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastp_annotations[row[0]] = row[12] + " E=" + str(row[10])

      ## Load pfam results
      print ("Loading pfam predictions from", input["pfam_results"], file=log_handle)
      with open(input["pfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 22)
          if (len(row) < 23): continue
          if row[3] not in pfam_annotations:
            pfam_annotations[row[3]] = row[1] + " " + row[22] + "E=" + str(row[6])

      ## Load rfam results
      print ("Loading rfam predictions from", input["rfam_results"], file=log_handle)
      with open(input["rfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 17)
          if (len(row) < 18): continue
          rfam_annotations[row[2]] = row[1] + " " + row[17] + "E=" + str(row[15])
  
      ## Load tmhmm results
      print ("Loading tmhmm predictions from", input["tmhmm_results"], file=log_handle)
      with open(input["tmhmm_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) > 2):
            tmhmm_annotations[row[0]] = row[1] + " " + row[2]
 
      ## Load deeploc results
      print ("Loading deeploc predictions from", input["deeploc_results"], file=log_handle)
      with open(input["deeploc_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 2): continue
          deeploc_annotations[row[0]] = str(row[1])

      ## Load targetp results
      translation =    {"SP": "Signal peptide",
                        "mTP": "Mitochondrial transit peptide",
                        "cTP": "chloroplast transit peptide",
                        "luTP": "thylakoidal lumen composite transit peptide",
                        "noTP": "no targeting peptide" }
      print("Loading targetp predictions from", input["targetp_results"], file=log_handle)
      with open(input["targetp_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = line.split()
          if (len(row) >= 8):
            targetp_annotations[row[0]] = translation[row[1]] + ", " + " ".join(row[7:])
          elif (len(row) >= 2):
            targetp_annotations[row[0]] = translation[row[1]]
      
      ## Do the work
      print ("Annotating FASTA file", input["transcriptome"], "to", output["transcriptome_annotated"], file=log_handle)
      with open(input["transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = record.id
          record.description = "TPM: " + expression_annotations.get(transcript_id)
          if transcript_id in blastx_annotations:
            record.description += "; blastx: " + blastx_annotations.get(transcript_id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")
      
      print ("Annotating FASTA file", input["proteome"], "to", output["proteome_annotated"], file=log_handle)
      with open(input["proteome"], "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = re.sub("\\.p[0-9]+$", "", record.id)
          record.description = "transdecoder: " + re.search("ORF type:([^,]+,score=[^,]+)", record.description).group(1)
          if transcript_id in expression_annotations:
            record.description += "; TPM: " + expression_annotations.get(transcript_id)
          if record.id in blastp_annotations:
            record.description += "; blastp: " + blastp_annotations.get(record.id)
          if record.id in pfam_annotations:
            record.description += "; pfam: " + pfam_annotations.get(record.id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          if record.id in tmhmm_annotations:
            record.description += "; tmhmm: " + tmhmm_annotations.get(record.id)
          if record.id in deeploc_annotations:
            record.description += "; deeploc: " + deeploc_annotations.get(record.id)
          if record.id in targetp_annotations:
            record.description += "; targetp: " + targetp_annotations.get(record.id)
          # Add sequence ID prefix from configuration
          if config["annotated_fasta_prefix"]:
            record.id = config["annotated_fasta_prefix"] + "|" + record.id
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")

      print ("Generating transcriptome_TPM_blast.csv table", file=log_handle)
      with open(input["expression"], "r") as input_csv_handle, open(output["tpm_blast_table"], "w") as output_csv_handle:
        csv_reader = csv.reader(input_csv_handle, delimiter="\t")
        csv_writer = csv.writer(output_csv_handle, delimiter=",")
        csv_columns = next(csv_reader)
        csv_columns[0] = "transcript"
        for i in range(1, len(csv_columns)):
          csv_columns[i] = "TPM(" + csv_columns[i] + ")"
        csv_columns.append("blastx")
        csv_writer.writerow(csv_columns)
        for row in csv_reader:
          row.append(blastx_annotations.get(row[0], ""))
          csv_writer.writerow(row)


