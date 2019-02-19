import os
import shutil
import re
import csv
import Bio.SeqIO
import Bio.Alphabet

shell.prefix("set -euo pipefail;")

configfile: "config.yaml"

TRINITY_HOME=os.path.dirname(shutil.which("Trinity"))

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis


rule all:
  input:
    config["samples_file"],
    "transcriptome.fasta",
    "transcriptome.pep",
    "transcriptome_stats.txt",
    "transcriptome_exN50.tsv.plot.pdf",
    "transcriptome_annotated.fasta",
    "transcriptome_annotated.pep",
    "transcriptome_TPM_blast.csv"


rule clean:
  shell:
    """
    rm -rf trinity_* tmp* log* TMHMM* kallisto* transcriptome* parallel* pipeliner* annotation*
    if [ -f {config[samples_file]} ]; then
      cut -f 2 < {config[samples_file]} | xargs --no-run-if-empty rm -rf
    fi
    """


rule trinity_inchworm_chrysalis:
  input:
    samples=config["samples_file"],
  output:
    "trinity_out_dir/recursive_trinity.cmds"
  log:
    "logs/log_trinity_inchworm_chrysalis.txt"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --samples_file {input} {config[trinity_parameters]} {config[strand_specific]} &> {log}
    """


checkpoint trinity_butterfly_split:
  input:
    "trinity_out_dir/recursive_trinity.cmds"
  output:
    directory("parallel/trinity_jobs")
  log:
    "logs/log_trinity_split.txt"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    mkdir -p {output} &> {log}
    split --numeric-suffixes=1 -l 100 {input} {output}/job_ &>> {log}
    """


rule trinity_butterfly_parallel:
  input:
    "parallel/trinity_jobs/job_{job_index}"
  output:
    "parallel/trinity_jobs/completed_{job_index}"
  log:
    "logs/log_trinity_parallel{job_index}.txt"
  params:
    memory="10"
  threads:
    1
  shell:
    """
    bash {input} &> {log}
    cp {input} {output} &>> {log}
    """

def trinity_completed_parallel_jobs(wildcards):
  parallel_dir = checkpoints.trinity_butterfly_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "job_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trinity_butterfly_merge:
  input:
    samples=config["samples_file"],
    jobs=trinity_completed_parallel_jobs
  output:
    cmds_completed="trinity_out_dir/recursive_trinity.cmds.completed",
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  log:
    "logs/log_trinity_butterfly_merge.txt"
  params:
    memory="200"
  threads:
    16
  shell:
    """
    cat {input.jobs} > {output.cmds_completed} 2> {log}
    Trinity --max_memory {params.memory}G --CPU {threads} --samples_file {input.samples} {config[trinity_parameters]} {config[strand_specific]} &>> {log}
    cp trinity_out_dir/Trinity.fasta {output.transcriptome} &>> {log}
    cp trinity_out_dir/Trinity.fasta.gene_trans_map {output.gene_trans_map} &>> {log}
    """


rule trinity_stats:
  input:
    transcriptome="transcriptome.fasta",
    expression="transcriptome_expression_isoform.tsv"
  output:
    stats="transcriptome_stats.txt",
    exN50="transcriptome_exN50.tsv",
    exN50plot="transcriptome_exN50.tsv.plot.pdf"
  log:
    "logs/log_trinity_exN50.txt"
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


rule transdecoder_longorfs:
  input:
    transcriptome="transcriptome.fasta",
  output:
    orfs="transcriptome.orfs"
  log:
    "logs/log_transdecoder_longorfs.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    rm -rf {input.transcriptome}.transdecoder_dir &> {log}
    TransDecoder.LongOrfs -t {input.transcriptome} --output_dir transdecoder &>> {log} 
    cp transdecoder/longest_orfs.pep {output.orfs} &>> {log}
    """


rule transdecoder_predict:
  input:
    transcriptome="transcriptome.fasta",
    pfam="annotations/pfamtransdecoder_orfs.out",
    blastp="annotations/sprotblastporfs_orfs.out"
  output:
    "transcriptome.pep"
  log:
    "logs/log_transdecoder_predict.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    TransDecoder.Predict -t {input.transcriptome} --output_dir transdecoder --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} &> {log}
    cp {input.transcriptome}.transdecoder.pep {output} &>> {log}
    """


rule align_reads:
  input:
    samples=config["samples_file"],
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  output:
    "RSEM_out.gene.TMM.EXPR.matrix"
  log:
    "logs/log_bowtie2.txt"
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
  input:
    samples=config["samples_file"],
    expression="kallisto.gene.counts.matrix"
  output:
    "edgeR_trans"
  log:
    "logs/log_trinity_DE.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    {TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input[1]} --samples_file {input[1]} --method edgeR --output {output} 2> {log}
    """


checkpoint fasta_split:
  input:
    "transcriptome.{extension}"
  output:
    directory("parallel/annotation_{extension}")
  log:
    "logs/log_{extension}_split.txt"
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
          filename = os.path.join(output[0], fileindex + "." + wildcards.extension)
          output_handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        if wildcards["extension"] == "pep":
          record.seq = record.seq.strip("*")
        # String the description, because some tools (e.g. deeploc) include it in their output
        record.description = ""
        Bio.SeqIO.write(record, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

def parallel_annotation_tasks(wildcards):
  parallel_dir = checkpoints.fasta_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}." + wildcards["extension"])).index
  completed_files = expand("parallel/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


rule annotation_merge_fasta:
  input:
    parallel_annotation_tasks
  output:
    "annotations/{task}_{extension}.out"
  log:
    "logs/log_{task}_{extension}_merge.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output} 2> {log}
    """


rule rfam_parallel:
  input:
    "parallel/annotation_fasta/{index}.fasta",
    "db/Rfam.cm"
  output:
    "parallel/rfam/{index}.out"
  log:
    "logs/log_rfam_{index}.txt"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    cmscan -E {config[e_value_threshold]} --rfam --cpu {threads} --tblout {output} {input[1]} {input[0]} &> {log}
    """

rule pfam_parallel:
  input:
    "parallel/annotation_pep/{index}.pep",
    "db/Pfam-A.hmm"
  output:
    "parallel/pfam/{index}.out"
  log:
    "logs/log_pfam_{index}.txt"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    hmmscan -E {config[e_value_threshold]} --cpu {threads} --tblout {output} {input[1]} {input[0]} &> {log}
    """

# Transdecoder requires --domtblout output
rule pfam_transdecoder_parallel:
  input:
    "parallel/annotation_orfs/{index}.orfs",
    "db/Pfam-A.hmm"
  output:
    "parallel/pfamtransdecoder/{index}.out"
  log:
    "logs/log_pfamtransdecoder_{index}.txt"
  params:
    memory="2"
  threads:
    2
  shell:
    """
    hmmscan -E {config[e_value_threshold]} --cpu {threads} --domtblout {output} {input[1]} {input[0]} &> {log}
    """


rule sprot_blastp_parallelorfs:
  input:
    "parallel/annotation_orfs/{index}.orfs",
    "db/uniprot_sprot.fasta"
  output:
    "parallel/sprotblastporfs/{index}.out"
  log:
    "logs/log_sprotblastporfs_{index}.txt"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastp -query {input[0]} -db {input[1]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule sprot_blastp_parallelpep:
  input:
    "parallel/annotation_pep/{index}.pep",
    "db/uniprot_sprot.fasta"
  output:
    "parallel/sprotblastppep/{index}.out"
  log:
    "logs/log_sprotblastppep_{index}.txt"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastp -query {input[0]} -db {input[1]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule sprot_blastx_parallel:
  input:
    "parallel/annotation_fasta/{index}.fasta",
    "db/uniprot_sprot.fasta"
  output:
    "parallel/sprotblastx/{index}.out"
  log:
    "logs/log_sprotblastx_{index}.txt"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastx -query {input[0]} -db {input[1]} -num_threads {threads} -evalue {config[e_value_threshold]} -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output} &> {log}
    """


rule tmhmm_parallel:
  input:
    "parallel/annotation_pep/{index}.pep"
  output:
    "parallel/tmhmm/{index}.out"
  log:
    "logs/log_tmhmm_{index}.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    tmhmm --short < {input} > {output} 2> {log}
    """


rule deeploc_parallel:
  input:
    "parallel/annotation_pep/{index}.pep"
  output:
    "parallel/deeploc/{index}.out"
  log:
    "logs/log_deeploc_{index}.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    deeploc -f {input} -o {output} &> {log}
    mv {output}.txt {output} &>> {log}
    """


rule kallisto:
  input:
    samples=config["samples_file"],
    transcriptome="transcriptome.fasta",
    gene_trans_map="transcriptome.gene_trans_map"
  output:
    "transcriptome_expression_isoform.tsv",
    "transcriptome_expression_gene.tsv"
  log:
    "logs/log_kallisto.txt"
  params:
    memory="2"
  threads:
    8
  shell:
    """
    {TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} {config[strand_specific]} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} &> {log}
    {TRINITY_HOME}/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input.gene_trans_map} */abundance.tsv &>> {log}
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp kallisto.isoform.TMM.EXPR.matrix {output[0]} &>> {log}
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp kallisto.isoform.TPM.not_cross_norm {output[0]} &>> {log}
    else
      echo Neither kallisto.isoform.TMM.EXPR.matrix or kallisto.isoform.TPM.not_cross_norm were produced
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp kallisto.gene.TMM.EXPR.matrix {output[1]} &>> {log}
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp kallisto.gene.TPM.not_cross_norm {output[1]} &>> {log}
    else
      echo Neither kallisto.gene.TMM.EXPR.matrix or kallisto.gene.TPM.not_cross_norm were produced
      exit 1
    fi
    """

 
rule download_sprot:
  output:
    "db/uniprot_sprot.fasta"
  log:
    "logs/log_download_sprot.txt"
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
  output:
    "db/Pfam-A.hmm"
  log:
    "logs/log_download_pfam.txt"
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
  output:
    "db/Rfam.cm"
  log:
    "logs/log_download_rfam.txt"
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
  output:
    "db/NOG.annotations.tsv"
  log:
    "logs/log_download_eggnog.txt"
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
  input:
    transcriptome="transcriptome.fasta",
    proteome="transcriptome.pep",
    expression="transcriptome_expression_isoform.tsv",
    blastx_results="annotations/sprotblastx_fasta.out",
    rfam_results="annotations/rfam_fasta.out",
    blastp_results="annotations/sprotblastppep_pep.out",
    pfam_results="annotations/pfam_pep.out",
    tmhmm_results="annotations/tmhmm_pep.out",
    deeploc_results="annotations/deeploc_pep.out"
  output:
    transcriptome_annotated="transcriptome_annotated.fasta",
    proteome_annotated="transcriptome_annotated.pep",
    tpm_blast_table="transcriptome_TPM_blast.csv"
  log:
    "logs/log_annotated_fasta.txt"
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
  
      ## Load kallisto results
      print ("Loading expression values from", input[2], file=log_handle)
      with open(input["expression"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter='\t')
        columns = next(csv_reader)
        for row in csv_reader:
          expression_annotations[row[0]] = columns[1] + "=" + str(row[1])
          for i in range(2, len(columns)):
            expression_annotations[row[0]] += " " + columns[i] + "=" + str(row[i])

      ## Load blastx results
      print ("Loading blastx results from", input[3], file=log_handle)
      with open(input["blastx_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastx_annotations[row[0]] = row[12] + " E=" + str(row[10])
  
      ## Load blastp results
      print ("Loading blastp results from", input[4], file=log_handle)
      with open(input["blastp_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 13): continue
          blastp_annotations[row[0]] = row[12] + " E=" + str(row[10])

      ## Load pfam results
      print ("Loading pfam predictions from", input[5], file=log_handle)
      with open(input["pfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 18)
          if (len(row) < 19): continue
          pfam_annotations[row[2]] = row[1] + " " + row[18] + "E=" + str(row[7])

      ## Load rfam results
      print ("Loading rfam predictions from", input[5], file=log_handle)
      with open(input["rfam_results"]) as input_handle:
        for line in input_handle:
          if (line.startswith("#")): continue
          row = re.split(" +", line, 17)
          if (len(row) < 18): continue
          rfam_annotations[row[2]] = row[1] + " " + row[17] + "E=" + str(row[15])
  
      ## Load tmhmm results
      print ("Loading tmhmm predictions from", input[6], file=log_handle)
      with open(input["tmhmm_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 6): continue
          tmhmm_annotations[row[0]] = row[2] + " " + row[3] + " " + row[4] + " " + row[5]
      
      ## Load deeploc results
      print ("Loading deeploc predictions from", input[7], file=log_handle)
      with open(input["deeploc_results"]) as input_handle:
        csv_reader = csv.reader(input_handle, delimiter="\t")
        for row in csv_reader:
          if (len(row) < 2): continue
          deeploc_annotations[row[0]] = str(row[1])
      
      ## Do the work
      print ("Annotating FASTA file", input[0], "to", output[0], file=log_handle)
      with open(input["transcriptome"], "r") as input_fasta_handle, open(output["transcriptome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = record.id
          record.description = "TPM: " + expression_annotations.get(transcript_id)
          if transcript_id in blastx_annotations:
            record.description += "; blastx: " + blastx_annotations.get(transcript_id)
          if transcript_id in rfam_annotations:
            record.description += "; rfam: " + rfam_annotations.get(transcript_id)
          Bio.SeqIO.write(record, output_fasta_handle, "fasta")
      
      print ("Annotating FASTA file", input[1], "to", output[1], file=log_handle)
      with open(input["proteome"], "r") as input_fasta_handle, open(output["proteome_annotated"], "w") as output_fasta_handle:
        for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
          transcript_id = re.sub("\\.p[0-9]+\$", "", record.id)
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


