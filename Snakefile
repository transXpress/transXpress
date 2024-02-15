import os
import sys
import shutil
import re
import csv
import Bio.SeqIO
import subprocess
from snakemake.utils import min_version
from Bio.Seq import Seq
import pandas as pd

min_version("5.4.1")

configfile: "config.yaml"

# for executable in ["samtools", "bowtie2", "kallisto", "signalp6", "targetp", "Trinity", "blastx", "blastp", "makeblastdb", "cmscan", "hmmscan", "fastqc", "rnaspades.py", "seqkit", "R"]:
#    if not shutil.which(executable):
#        sys.stderr.write("Warning: Cannot find %s in your PATH\n" % (executable))

#TRINITY_EXECUTABLE_PATH=shutil.which("Trinity")
#TRINITY_HOME=os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))) ##Have to resolve the symbolic link that conda makes

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis

# These rules don't need to be sent to a cluster
localrules: all, clean, trimmomatic_split, trimmomatic_merge, samples_yaml_conversion, trinity_butterfly_split, transcriptome_copy, compare_qc_after_trim, fastqc_before_trim_warnings, IGV


rule all:
  """
  List of target files of the transxpress pipeline.
  """
  input: 
    "samples_trimmed.txt",
    "transcriptome.fasta",
    "transcriptome.pep",
    "transcriptome_clst.fasta",
    "transcriptome_stats.txt",
    "ExN50_plot.pdf",
    "transcriptome_annotated.fasta",
    "transcriptome_annotated.pep",
    "transcriptome_TPM_blast.csv",
    "busco",
    "busco_report.txt",
    "fastqc_before_trim",
    "multiqc_before_trim",
    "fastqc_after_trim",
    "multiqc_after_trim",
    "multiqc_before_trim.txt",
    "multiqc_after_trim.txt",
    "WARNING_fastqc_before_trim_overview.txt",
    "FastQC_comparison_after_trim.txt",
    "edgeR_trans",
    ##"clusters.tsv" - can only uncomment if CD-HIT option is set to true in the config.yaml file

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


rule fastqc_before_trim:
  """
  Runs fastQC on individual input reads files.
  """
  input:
    samples=config["samples_file"]
  output:
    directory("fastqc_before_trim")
  log:
    "logs/fastqc_before_trim.log"
  conda:
    "envs/qc.yaml"
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


rule multiqc_before_trim:
  """
  Creates multiQC report from individual fastQC reports.
  """
  input:
    "fastqc_before_trim"
  output:
    out_dir=directory("multiqc_before_trim"),
    report="multiqc_before_trim.txt"
  log:
    "logs/multiqc_before_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_dir} &> {log}
    multiqc -o {output.out_dir} {input} &>> {log}
    cp multiqc_before_trim/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule fastqc_before_trim_warnings:
  """
  Prints warnings for FastQC modules which produced some warnings or fails of the input data.
  """
  input:
    fastqc_report="multiqc_before_trim.txt"
  output:
    warning_file="WARNING_fastqc_before_trim_overview.txt"
  log:
    "logs/fastqc_before_trim_warnings.log"
  run:

    mapping = {
    "basic_statistics": "BASIC STATISTICS",
    "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
    "per_tile_sequence_quality": "PER TILE SEQUENCE QUALITY",
    "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
    "per_base_sequence_content": "PER BASE SEQUENCE CONTENT",
    "per_sequence_gc_content": "PER SEQUENCE GC CONTENT",
    "per_base_n_content": "PER BASE N CONTENT",
    "sequence_length_distribution": "SEQUENCE LENGTH DISTRIBUTION",
    "sequence_duplication_levels": "SEQUENCE DUPLICATION LEVELS",
    "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
    "adapter_content": "ADAPTER CONTENT"
    }

    reasons_mapping = {
    "basic_statistics": "BASIC STATISTICS",
    
    "per_base_sequence_quality":
    """A warning is issued if the lower quartile for any base is less than 10, or if the median for any base is less than 25.
This module will raise a failure if the lower quartile for any base is less than 5 or if the median for any base is less than 20.

Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",
    
    "per_tile_sequence_quality":
"""This module issues a warning if any tile shows a mean Phred score more than 2 less than the mean for that base across all tiles.
This module will issue a warning if any tile shows a mean Phred score more than 5 less than the mean for that base across all tiles.

Common reasons for warnings:
Whilst warnings in this module can be triggered by individual specific events we have also observed that greater variation in the phred scores attributed to tiles can also appear when a flowcell is generally overloaded. In this case events appear all over the flowcell rather than being confined to a specific area or range of cycles. We would generally ignore errors which mildly affected a small number of tiles for only 1 or 2 cycles, but would pursue larger effects which showed high deviation in scores, or which persisted for several cycles.
""",
    
    "per_sequence_quality_scores": 
"""A warning is raised if the most frequently observed mean quality is below 27 - this equates to a 0.2% error rate.
An error is raised if the most frequently observed mean quality is below 20 - this equates to a 1% error rate.

Common reasons for warnings:
This module is generally fairly robust and errors here usually indicate a general loss of quality within a run. For long runs this may be alleviated through quality trimming. If a bi-modal, or complex distribution is seen then the results should be evaluated in concert with the per-tile qualities (if available) since this might indicate the reason for the loss in quality of a subset of sequences.
""",
    
    "per_base_sequence_content": 
"""Common reasons for warnings:
The most common reason for warnings and failures in this module is a general degradation of quality over the duration of long runs. In general sequencing chemistry degrades with increasing read length and for long runs you may find that the general quality of the run falls to a level where a warning or error is triggered.

If the quality of the library falls to a low level then the most common remedy is to perform quality trimming where reads are truncated based on their average quality. For most libraries where this type of degradation has occurred you will often be simultaneously running into the issue of adapter read-through so a combined adapter and quality trimming step is often employed.

Another possibility is that a warn / error is triggered because of a short loss of quality earlier in the run, which then recovers to produce later good quality sequence. This can happen if there is a transient problem with the run (bubbles passing through a flowcell for example). You can normally see this type of error by looking at the per-tile quality plot (if available for your platform). In these cases trimming is not advisable as it will remove later good sequence, but you might want to consider masking bases during subsequent mapping or assembly.

If your library has reads of varying length then you can find a warning or error is triggered from this module because of very low coverage for a given base range. Before committing to any action, check how many sequences were responsible for triggering an error by looking at the sequence length distribution module results.
""",
    
    "per_sequence_gc_content": 
"""A warning is raised if the sum of the deviations from the normal distribution represents more than 15% of the reads.
This module will indicate a failure if the sum of the deviations from the normal distribution represents more than 30% of the reads.

Common reasons for warnings:
Warnings in this module usually indicate a problem with the library. Sharp peaks on an otherwise smooth distribution are normally the result of a specific contaminant (adapter dimers for example), which may well be picked up by the overrepresented sequences module. Broader peaks may represent contamination with a different species.
""",
    
    "per_base_n_content": 
"""This module raises a warning if any position shows an N content of >5%.
This module will raise an error if any position shows an N content of >20%.

Common reasons for warnings:
The most common reason for the inclusion of significant proportions of Ns is a general loss of quality, so the results of this module should be evaluated in concert with those of the various quality modules. You should check the coverage of a specific bin, since it is possible that the last bin in this analysis could contain very few sequences, and an error could be prematurely triggered in this case.

Another common scenario is the incidence of a high proportions of N at a small number of positions early in the library, against a background of generally good quality. Such deviations can occur when you have very biased sequence composition in the library to the point that base callers can become confused and make poor calls. This type of problem will be apparent when looking at the per-base sequence content results.
""",
    
    "sequence_length_distribution": 
"""This module will raise a warning if all sequences are not the same length.
This module will raise an error if any of the sequences have zero length.

Common reasons for warnings:
For some sequencing platforms it is entirely normal to have different read lengths so warnings here can be ignored.""",
    
    "sequence_duplication_levels": 
"""This module will issue a warning if non-unique sequences make up more than 20% of the total.
This module will issue a error if non-unique sequences make up more than 50% of the total.

Common reasons for warnings:
The underlying assumption of this module is of a diverse unenriched library. Any deviation from this assumption will naturally generate duplicates and can lead to warnings or errors from this module.
        
In general there are two potential types of duplicate in a library, technical duplicates arising from PCR artefacts, or biological duplicates which are natural collisions where different copies of exactly the same sequence are randomly selected. From a sequence level there is no way to distinguish between these two types and both will be reported as duplicates here.
        
A warning or error in this module is simply a statement that you have exhausted the diversity in at least part of your library and are re-sequencing the same sequences. In a supposedly diverse library this would suggest that the diversity has been partially or completely exhausted and that you are therefore wasting sequencing capacity. However in some library types you will naturally tend to over-sequence parts of the library and therefore generate duplication and will therefore expect to see warnings or error from this module.
        
In RNA-Seq libraries sequences from different transcripts will be present at wildly different levels in the starting population. In order to be able to observe lowly expressed transcripts it is therefore common to greatly over-sequence high expressed transcripts, and this will potentially create large set of duplicates. This will result in high overall duplication in this test, and will often produce peaks in the higher duplication bins. This duplication will come from physically connected regions, and an examination of the distribution of duplicates in a specific genomic region will allow the distinction between over-sequencing and general technical duplication, but these distinctions are not possible from raw fastq files. A similar situation can arise in highly enriched ChIP-Seq libraries although the duplication there is less pronounced. Finally, if you have a library where the sequence start points are constrained (a library constructed around restriction sites for example, or an unfragmented small RNA library) then the constrained start sites will generate huge dupliction levels which should not be treated as a problem, nor removed by deduplication. In these types of library you should consider using a system such as random barcoding to allow the distinction of technical and biological duplicates.
""",
    
    "overrepresented_sequences": 
"""This module will issue a warning if any sequence is found to represent more than 0.1% of the total.
This module will issue an error if any sequence is found to represent more than 1% of the total.
        
Common reasons for warnings:
This module will often be triggered when used to analyse small RNA libraries where sequences are not subjected to random fragmentation, and the same sequence may natrually be present in a significant proportion of the library.
""",
    
    "adapter_content": 
"""This module will issue a warning if any sequence is present in more than 5% of all reads.
This module will issue an error if any sequence is present in more than 10% of all reads.
        
Common reasons for warnings:
Any library where a reasonable proportion of the insert sizes are shorter than the read length will trigger this module. This does not indicate a problem as such - just that the sequences will need to be adapter trimmed before proceeding with any downstream analysis.
"""
    }



    with open(output["warning_file"], "w") as out_file:
      df1 = pd.read_csv(input["fastqc_report"], sep='\t')
      df1 = df1.set_index('Sample')

      # filter only columns containing 'pass','warn' or 'fail'
      df1_filtered = df1.loc[:, df1.isin(['pass','warn','fail']).all()]
      df1_filtered = df1_filtered.loc[:, df1_filtered.columns != 'basic_statistics'] # basic statistics module always passes

      n_samples = df1.shape[0]

      for column in df1_filtered.columns:
        if column in mapping.keys():
          out_file.writelines(mapping[column]+'\n')
        else:
          out_file.writelines(column+'\n')
        results = df1_filtered[column].tolist()
        if ('warn' in results) or ('fail' in results):
          n_warn = results.count('warn')
          n_fail = results.count('fail')
          n_warn_perc = round(n_warn/n_samples*100,2)
          n_fail_perc = round(n_fail/n_samples*100,2)
          if n_warn > 0:
            out_file.writelines(f'{n_warn} out of {n_samples} produced "WARN" ({n_warn_perc}%)\n')
          if n_fail > 0:
            out_file.writelines(f'{n_fail} out of {n_samples} produced "FAIL" ({n_fail_perc}%)\n')
          if column in reasons_mapping.keys():
            out_file.writelines('\n' + reasons_mapping[column])
        else:
          out_file.writelines('ok'+'\n')
      out_file.writelines('source: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/\n')
    
    # print to stdout
    with open(output["warning_file"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in: WARNING_fastqc_before_trim_overview.txt')
    
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
  conda:
    "envs/trimmomatic.yaml"
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

rule fastqc_after_trim:
  """
  Runs fastQC on individual trimmed reads files.
  """
  input:
    samples="samples_trimmed.txt"
  output:
    directory("fastqc_after_trim")
  log:
    "logs/fastqc_after_trim.log"
  conda:
    "envs/qc.yaml"
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


rule multiqc_after_trim:
  """
  Creates multiQC report from individual fastQC reports after the trimming.
  """
  input:
    "fastqc_after_trim"
  output:
    out_directory=directory("multiqc_after_trim"),
    report="multiqc_after_trim.txt"
  log:
    "logs/multiqc_after_trim.log"
  conda:
    "envs/qc.yaml"
  params:
    memory="4"
  threads:
    1
  shell:
    """
    # common error resolved by those two export commands
    export LC_ALL=C.UTF-8
    export LANG=C.UTF-8
    mkdir {output.out_directory} &> {log}
    multiqc -o {output.out_directory} {input} &>> {log}
    cp multiqc_after_trim/multiqc_data/multiqc_fastqc.txt {output.report} &>> {log}
    """

rule compare_qc_after_trim:
  """
  Compares the quality of reads before and after trimming.
  Compares PER BASE SEQUENCE QUALITY, PER SEQUENCE QUALITY SCORES, 
  OVERREPRESENTED SEQUENCES and ADAPTER CONTENT modules of FastQC.
  """
  input:
    before="multiqc_before_trim.txt",
    after="multiqc_after_trim.txt"
  output:
    result="FastQC_comparison_after_trim.txt"
  log:
    "logs/compare_qc_after_trim.log"
  run:

    mapping = {
      "per_base_sequence_quality": "PER BASE SEQUENCE QUALITY",
      "per_sequence_quality_scores": "PER SEQUENCE QUALITY SCORES",
      "overrepresented_sequences": "OVERREPRESENTED SEQUENCES",
      "adapter_content": "ADAPTER CONTENT"
    }

    def compare_module(module_name, before_trim, after_trim, out_file):
      module_name = mapping[module_name]
      out_file.writelines(f"""**{module_name}**
Before trimming:
pass: {str(before_trim.count('pass'))} out of {str(len(before_trim))} ({str(before_trim.count('pass')/len(before_trim)*100)}%)
warn: {str(before_trim.count('warn'))} out of {str(len(before_trim))} ({str(before_trim.count('warn')/len(before_trim)*100)}%)
fail: {str(before_trim.count('fail'))} out of {str(len(before_trim))} ({str(before_trim.count('fail')/len(before_trim)*100)}%)
      
After trimming:
pass: {str(after_trim.count('pass'))} out of {str(len(after_trim))} ({str(after_trim.count('pass')/len(after_trim)*100)}%)
warn: {str(after_trim.count('warn'))} out of {str(len(after_trim))} ({str(after_trim.count('warn')/len(after_trim)*100)}%)
fail: {str(after_trim.count('fail'))} out of {str(len(after_trim))} ({str(after_trim.count('fail')/len(after_trim)*100)}%)
{'-' * 50}
""")

    before = pd.read_csv(input["before"], sep='\t')
    before = before.set_index('Sample')

    after = pd.read_csv(input["after"], sep='\t')
    after = after.set_index('Sample')

    modules = ['per_base_sequence_quality', 'per_sequence_quality_scores', 'overrepresented_sequences', 'adapter_content']

    with open(output["result"], "w") as out_file:
      out_file.write("FASTQC COMPARISON BEFORE AND AFTER TRIMMING\n")
      out_file.write('-' * 50 + '\n')

      for module in modules:
        before_list = before[module].tolist()
        after_list = after[module].tolist()
        compare_module(module, before_list, after_list, out_file)

    # print the results also in the terminal
    with open(output["result"], "r") as out_file:
      for line in out_file.readlines():
        print(line.strip())
      print('See this output in FastQC_comparison_after_trim.txt')

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
  conda:
    "envs/trinity_utils.yaml"
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
  conda:
    "envs/trinity_utils.yaml"
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
    transcriptome="trinity_out_dir.Trinity.fasta", #new output files live in base directory now, not trinity_out_dir.
    gene_trans_map="trinity_out_dir.Trinity.fasta.gene_trans_map" #according to trinity dev, this is a feature, not a bug
  log:
    "logs/trinity_final.log"
  conda:
    "envs/trinity_utils.yaml"
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
  conda:
    "envs/rnaspades.yaml"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    rnaspades.py --dataset {input.samples} -t {threads} -m {params.memory} -o rnaspades_out {config[strand_specific]} --only-assembler {config[kmers]} &> {log}
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
    gene_trans_map="transcriptome.gene_trans_map",
    #create copies the same way older Trinity versions did, for parity, in case other dependencies on that path exist?
    redundant_transcriptome="trinity_out_dir/Trinity.fasta",
    redundant_gene_trans_map="trinity_out_dir/Trinity.fasta.gene_trans_map"
  log:
    "logs/transcriptome_copy.log"
  shell:
    """
    cp -p {input.transcriptome} {output.transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.gene_trans_map} &>> {log}
    #create copies the same way older Trinity versions did, for parity, in case other dependencies on that path exist?
    cp -p {input.transcriptome} {output.redundant_transcriptome} &> {log}
    cp -p {input.gene_trans_map} {output.redundant_gene_trans_map} &>> {log}
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
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    assembler={config[assembler]}
    if [ "$assembler" = 'trinity' ]; then
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} > {output.stats} 2> {log}
    else
        $TRINITY_HOME/util/TrinityStats.pl {input.transcriptome} | sed -e 's/trinity/rnaspades/g' > {output.stats} 2> {log}
    fi
    $TRINITY_HOME/util/misc/contig_ExN50_statistic.pl {input.expression} {input.transcriptome} > {output.exN50} 2>> {log}
    $TRINITY_HOME/util/misc/plot_ExN50_statistic.Rscript {output.exN50} &>> {log}
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
  conda:
    "envs/busco.yaml"
  params:
    memory="10"
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

    status=$?

    if [ $status -eq 0 ]
    then
      echo "BUSCO run completed successfully" &>> {log}
      cp busco/short_summary*.txt {output.report} &>> {log}
    else
      echo "BUSCO run failed" &>> {log}
      exit 1 &>> {log}
    fi
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
  conda:
    "envs/transdecoder.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    rm -rf {input.transcriptome}.transdecoder_dir &> {log}
    TransDecoder.LongOrfs -t {input.transcriptome} --output_dir transdecoder &>> {log} 
    cp -p transdecoder/{input.transcriptome}.transdecoder_dir/longest_orfs.pep {output.orfs} &>> {log}
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
  conda:
    "envs/transdecoder.yaml"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    TransDecoder.Predict -t {input.transcriptome} --output_dir transdecoder --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} &> {log}
    cp -p transdecoder/{input.transcriptome}.transdecoder.pep {output} &>> {log}
    """

rule cd_hit:
  """
  (optional)
  Apply CD-HIT on transcriptome.pep to get rid of redundant sequences (saves space and time complexity). CD-HIT reduces space by clustering 100% identical 
  sequences into one representative sequence. If config[cd-hit] is set to False, we just copy the input file into the output file.
  output are two files: a fasta file of representative sequences and a text file of list of clusters
  """
  input:
    "transcriptome.pep"
  output:
    "transcriptome_clst.pep"
  log:
    "logs/cd_hit.log"
  conda:
    "envs/cdhit.yaml"
  params:
    memory="8"
  threads:
    1
  shell:
    """
    if [ {config[cd-hit]} = 'true' ]; then
      cd-hit -i {input} -o {output} -c 1.00 -n 5 -g 1 -d 0 &>> {log}
    else
      cp {input} {output}
    fi
    """
  
rule parse_cdhit:
  """
  Parse the CD-HIT output.
  """
  input:
    clustered="transcriptome_clst.pep.clstr" # if cd-hit is not run by the pipeline this file will be missing, but this rule can only be run if the user wishes to have the tsv file in the outputs of the pipeline
  output:
    tsv="clusters.tsv"
  log:
    "logs/parse_cdhit.log"
  threads:
    1
  run:
    cl_size = 0          
    clstr_cl_size = ""  
    this_cluster = []    

    def process_this_cluster(this_cluster, output):
      t = sorted(this_cluster, key=lambda x: (x[2], x[1]), reverse=True)
      longest = 0
      for i in t:
        if i[2]:
          longest = i[1]
      for i in t:
        if i[1] is None:
          continue
        else:
          cov = int(i[1] / longest * 100) 
        print(f"{i[0]}\t{clstr_cl_size}\t{cl_size}\t{i[1]}\t{i[2]}\t{i[3]}\t{cov}%", file=output)

    with open(output["tsv"], "w") as output_file:
      print("id\tclstr\tclstr_size\tlength\tclstr_rep\tclstr_iden\tclstr_cov", file=output_file)  
      for line in open(input["clustered"]):
        if line.startswith('>'):
          if cl_size > 0:
            process_this_cluster(this_cluster, output_file) 
          if line.startswith('>Cluster '):
            clstr_cl_size = line.split()[1]  
          cl_size = 0
          this_cluster = []
        else:
          id, len, rep, iden = None, None, None, None
          match = re.match(r'\d+\t(\d+)[a-z]{2}, >NODE_\d+_length_\d+_cov_\d+\.\d+_(g\d+_i\d+\.p\d+)', line)
          if match:
            len = int(match.group(1))
            id = match.group(2)
            rep = 1   
            iden = 100 
          elif re.match(r'\d+\t(\d+)[a-z]{2}, >(.+)\.\.\.', line):
            match = re.match(r'\d+\t(\d+)[a-z]{2}, >(.+)\.\.\.', line)
            len = int(match.group(1))
            id = match.group(2)
            rep = 0  
            match = re.search(r'(\d+%|\d+\.\d+%)$', line)
            iden = match.group(1)
          else:
            print("***********")
          this_cluster.append((id, len, rep, iden))
          cl_size += 1

      if cl_size > 0:
        process_this_cluster(this_cluster, output_file)


rule reduce_transcriptome: 
  """
  Combines output of CD-HIT and the transcriptome.fasta file to a clustered transcriptome_clst.fasta file. Output will be a reduced 
  transcriptome which will be used as the input for kallisto.
  """
  input:
    clustered_proteome="transcriptome_clst.pep",
    transcriptome_path="transcriptome.fasta"
  output:
    transcriptome_reduced="transcriptome_clst.fasta"
  log:
    "logs/transcriptome_clust.log"
  threads:
    1
  params:
    memory="8"
  run:
    cd_hit_performed = config["cd-hit"]
    if cd_hit_performed == "true":
    	cd_hit_performed = True
    else:
	    cd_hit_performed = False
	
    if cd_hit_performed:
      headers = []

      with open(input["clustered_proteome"], 'r') as proteome:
        for record in Bio.SeqIO.parse(proteome, "fasta"):
          headers.append(record.id)

      headers_short = [header.split('_') for header in headers]
      headers_short = ['_'.join(sublist[6:]) for sublist in headers_short]
      headers_short = [header_short.split('.')[0] for header_short in headers_short]

      with open(input["transcriptome_path"], 'r') as transcriptome, open(output["transcriptome_reduced"], 'w') as transcriptome_reduced:
        for record in Bio.SeqIO.parse(transcriptome, 'fasta'):
          id_short = record.id.split('_')
          id_short = '_'.join(id_short[6:])
        if id_short in headers_short:
          Bio.SeqIO.write(record, transcriptome_reduced, 'fasta-2line')
    else:
      shutil.copy(input["transcriptome_path"], output["transcriptome_reduced"])


checkpoint align_reads:
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
    "transcriptome.fasta.bowtie2.ok"
  log:
    "logs/bowtie2.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="100"
  threads:
    16
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')


    assembler="{config[assembler]}"
    strand_specific="{config[strand_specific]}"
    
    if [ $assembler = "rnaspades" ]
    then
      if [[ $strand_specific = "--ss rf" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input[1]} --SS_lib_type RF --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map {input[2]} &> {log}
      elif [[ $strand_specific = "--ss fr" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input[1]} --SS_lib_type FR --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map {input[2]} &>> {log}
      else
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input[1]} --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map {input[2]} &>> {log}
      fi
    else
      $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts {input[1]} {config[strand_specific]} --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map {input[2]} &>> {log}
    fi
    """

checkpoint prepare_samples_for_IGV:
  """
  Not run by default.
  Moves bam alignment files to bowtie_alignments directory
  and prepares the files to be used by IGV by sorting and 
  indexing them.
  """
  input:
    alignment="{sample}/bowtie2.bam"
  output:
    alignment="bowtie_alignments/{sample}.bam",
    sorted_alignment="bowtie_alignments/{sample}.sorted.bam",
    indexed_alignment="bowtie_alignments/{sample}.sorted.bam.bai"
  log:
    "logs/prepare_for_IGV_{sample}.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cp {input.alignment} {output.alignment} &> {log}
    samtools sort --threads {threads} -o {output.sorted_alignment} {output.alignment} &>> {log}
    samtools index {output.sorted_alignment} &>> {log}
    """ 

def get_igv(wildcards):
  indeces = glob_wildcards("{sample}/bowtie2.bam").sample
  completed = expand(os.path.join("bowtie_alignments", "{sample}.sorted.bam.bai"),sample=indeces)
  return completed

rule IGV:
  """
  Not run by default.
  Rule to be called to get all bowtie alignments ready for IGV.
  """
  input:
    get_igv,
    "transcriptome.fasta.bowtie2.ok"
  output:
    "igv.ok"
  log:
    "logs/IGV_combine.log"
  shell:
    """
    touch {output}
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
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')


    num_replicates=`awk '{{print $2}}' {input.samples} | sort | uniq | wc -l` &> {log}
    num_samples=`awk '{{print $1}}' {input.samples} | sort | uniq | wc -l` &>> {log}
    num_replicates_minus_samples=$((num_replicates - num_samples)) &>> {log}
    if [ $num_replicates_minus_samples -gt 1 ] &>> {log}
    then &>> {log}
	    $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} &>> {log}
    else &>> {log}
	    echo "No biological replicates to run proper differential expression analysis, last-resorting to edgeR with --dispersion 0.1" &>> {log}
      $TRINITY_HOME/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input.expression} --method edgeR --output {output} --samples_file {input.samples} --dispersion {config[dispersion]} &>> {log}
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
    "transcriptome_clst.pep"
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
  conda:
    "envs/rfam.yaml"
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
  conda:
    "envs/pfam.yaml"
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
  conda:
    "envs/blast.yaml"
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
  conda:
    "envs/blast.yaml"
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

rule signalp_parallel:
  """
  Runs signalp on smaller protein files in parallel to predict signal peptides.
  """
  input:
    "annotations/chunks_pep/{index}.pep"
  output:
    "annotations/signalp/{index}.out"
  log:
    "logs/signalp_{index}.log"
  params:
    memory="8"
  threads:
    1
  shell:
    """
    mkdir -p annotations/signalp &> {log}
    signalp6 --fastafile {input} --organism eukarya --output_dir signalp_{wildcards.index} --format none --mode fast &>> {log}
    mv signalp_{wildcards.index}/prediction_results.txt {output}
    rm -r signalp_{wildcards.index}
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
    gene_trans_map="transcriptome.gene_trans_map",
    clustered_transcriptome="transcriptome_clst.fasta"
  output:
    "transcriptome_expression_isoform.tsv",
    "transcriptome_expression_gene.tsv",
    "kallisto.gene.counts.matrix"
  log:
    "logs/kallisto.log"
  conda:
    "envs/trinity_utils.yaml"
  params:
    memory="16" # increased memory from 2 to 8 since it was not sufficient
  threads:
    8
  shell:
    """
    TRINITY_HOME=$(python -c 'import os;import shutil;TRINITY_EXECUTABLE_PATH=shutil.which("Trinity");print(os.path.dirname(os.path.join(os.path.dirname(TRINITY_EXECUTABLE_PATH), os.readlink(TRINITY_EXECUTABLE_PATH))))')

    assembler="{config[assembler]}"
    strand_specific="{config[strand_specific]}"
	
    if [ {config[cd-hit]} = "true" ]; then
        transcriptome={input.clustered_transcriptome}
    else
        transcriptome={input.transcriptome}
    fi

    if [ $assembler = "rnaspades" ]
    then
      if [[ $strand_specific = "--ss rf" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $transcriptome --SS_lib_type RF --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
      elif [[ $strand_specific = "--ss fr" ]]
      then
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $transcriptome --SS_lib_type FR --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
      else
        $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $transcriptome --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
      fi
    else
      $TRINITY_HOME/util/align_and_estimate_abundance.pl --transcripts $transcriptome {config[strand_specific]} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
    fi
    
    $TRINITY_HOME/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input.gene_trans_map} */abundance.tsv &>> {log}
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
  conda:
    "envs/blast.yaml"
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
  conda:
    "envs/pfam.yaml"
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
  conda:
    "envs/rfam.yaml"
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
    clustered_transcriptome="transcriptome_clst.fasta",
    proteome="transcriptome.pep",
    clustered_proteome="transcriptome_clst.pep",
    expression="transcriptome_expression_isoform.tsv",
    blastx_results="annotations/sprotblastx_fasta.out",
    rfam_results="annotations/rfam_fasta.out",
    blastp_results="annotations/sprotblastp_orfs.out",
    pfam_results="annotations/pfam_orfs.out",
    tmhmm_results="annotations/tmhmm_pep.out",
    signalp_results="annotations/signalp_pep.out",
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
      signalp_annotations = {}
      targetp_annotations = {}

      cd_hit_performed = config["cd-hit"]
      if cd_hit_performed == "true":
        cd_hit_performed = True
      else:
        cd_hit_performed = False
  
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

      ## Load signalp results
      print ("Loading signalp predictions from", input["signalp_results"], file=log_handle)
      with open(input["signalp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 2): continue
          prediction = str(row[1])
          if prediction == "OTHER":
            signalp_annotations[row[0]] = "OTHER (No SP)"
          elif prediction == "SP":
            signalp_annotations[row[0]] = "standard secretory signal peptides (Sec/SPI)"
          elif prediction == "LIPO":
            signalp_annotations[row[0]] = "lipoprotein signal peptides (Sec/SPII)"
          elif prediction == "TAT":
            signalp_annotations[row[0]] = "Tat signal peptides (Tat/SPI)"
          elif prediction == "TATLIPO":
            signalp_annotations[row[0]] = "Tat lipoprotein signal peptides (Tat/SPII)"
          elif prediction == "PILIN":
            signalp_annotations[row[0]] = "Pilin and pilin-like signal peptides (Sec/SPIII)"
          else:
            signalp_annotations[row[0]] = str(row[1])

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
      if cd_hit_performed:
        input_transcriptome = input["clustered_transcriptome"]
      else:
        input_transcriptome = input["transcriptome"]
      with open(input["input_transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = record.id
          record.description = "TPM: " + expression_annotations.get(transcript_id)
          if transcript_id in blastx_annotations:
            record.description += "; blastx: " + blastx_annotations.get(transcript_id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")
      
      print ("Annotating FASTA file", input["proteome"], "to", output["proteome_annotated"], file=log_handle)
      if cd_hit_performed:
        input_proteome = input["clustered_proteome"]
      else:
        input_proteome = input["proteome"]
      with open(input_proteome, "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
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
          if record.id in signalp_annotations:
            record.description += "; signalp: " + signalp_annotations.get(record.id)
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


