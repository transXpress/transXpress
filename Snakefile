# Note: must avoid spaces in filenames in samples.txt

import os
import shutil
import re
import csv
import Bio.SeqIO
import Bio.Alphabet

shell.prefix("set -euo pipefail;")

configfile: "config.yaml"


# Note: there is another signalp in /usr/local/bin, but that one will crash with > 10000 sequences
SIGNALP_ORGANISM="euk"


TRINITY_HOME=os.path.dirname(shutil.which("Trinity"))

STRAND_SPECIFIC = "" # --SS_lib_type RF

TRINITY_PARAMS = " --seqType fq"
TRINITY_PARAMS += " --trimmomatic --quality_trimming_params \"ILLUMINACLIP:" + TRINITY_HOME + "/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:3:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:25\""
TRINITY_PARAMS += " " + STRAND_SPECIFIC
 
# may consider changing TRINITY_PARAMS += "  --min_kmer_cov 1 --min_glue 2  --no_normalize_reads"

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
    rm -rf trinity_* tmp* log* TMHMM* tmhmm* kallisto* transcriptome* proteome* pfam* sprot* signalp* parallel* pipeliner* annotation*
    if [ -f samples.txt ]; then
      cut -f 2 < samples.txt | xargs --no-run-if-empty rm -rf
    fi
    """


rule trinity_inchworm_chrysalis:
  input:
    "samples.txt"
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
    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --samples_file {input} {TRINITY_PARAMS} 2>&1 > {log}
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
    mkdir -p {output} 2> {log}
    split -l 100 {input} {output}/job_ 2> {log}
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
    bash {input} 2> {log}
    cp {input} {output} 2> {log}
    """

def trinity_completed_parallel_jobs(wildcards):
  parallel_dir = checkpoints.trinity_butterfly_split.get(**wildcards).output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "job_{job_index}")).job_index
  completed_ids = expand(os.path.join(parallel_dir,"completed_{job_index}"), job_index=job_ids)
  return completed_ids

rule trinity_butterfly_merge:
  input:
    samples="samples.txt",
    cmds="trinity_out_dir/recursive_trinity.cmds",
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
    cp -n {input.cmds} {output.cmds_completed} 2> {log}
    Trinity --max_memory {params.memory}G --CPU {threads} --samples_file {input.samples} {TRINITY_PARAMS} 2>&1 >> {log}
    mv trinity_out_dir/Trinity.fasta {output.transcriptome} 2>> {log}
    mv trinity_out_dir/Trinity.fasta.gene_trans_map {output.gene_trans_map} 2>> {log}
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
    {TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript {output.exN50} 2>> {log}
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
    rm -rf {input.transcriptome}.transdecoder_dir
    TransDecoder.LongOrfs -t {input.transcriptome} 2>&1 > {log} 
    cp transcriptome.fasta.transdecoder_dir/longest_orfs.pep {output.orfs} 2>&1 >> {log}
    """


rule transdecoder_predict:
  input:
    transcriptome="transcriptome.fasta",
    pfam="annotations_orfs/pfam.out",
    blastp="annotations_orfs/sprot_blastp.out"
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
    TransDecoder.Predict -t {input.transcriptome} --retain_pfam_hits {input.pfam} --retain_blastp_hits {input.blastp} 2>&1 > {log}
    mv {input.transcriptome}.transdecoder.cds {output} 2>&1 >> {log}
    """


rule align_reads:
  input:
    "samples.txt",
    "transcriptome.fasta",
    "transcriptome.gene_trans_map"
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
    {TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts {input[1]} {STRAND_SPECIFIC} --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method RSEM --aln_method bowtie2 --gene_trans_map {input[2]}
    # samtools sort -l 9 -@ {threads} -T  ??? ???/bowtie2.bam > bowtie2.sorted.bam
    # samtools index bowtie2.sorted.bam
    """


rule trinity_DE:
  input:
    "samples.txt",
    "kallisto.gene.counts.matrix"
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
      parser = Bio.SeqIO.parse(input_handle, "fasta")
      output_handle = None
      count = 0
      for entry in parser:
        if (count % 1000 == 0):
          if output_handle is not None:
            output_handle.close()
          fileindex = str(int(count / 1000) + 1);
          filename = os.path.join(output[0], fileindex + "." + wildcards.extension)
          print("Writing to file " + filename)
          output_handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        entry.seq = entry.seq.strip("*")
        Bio.SeqIO.write(entry, output_handle, "fasta")
        count += 1
      if output_handle is not None:
        output_handle.close()
        

def parallel_annotation_orfs(wildcards):
  parallel_dir = checkpoints.fasta_split.get(extension="orfs").output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.orfs")).index
  completed_files = expand("parallel/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files

def parallel_annotation_fasta(wildcards):
  parallel_dir = checkpoints.fasta_split.get(extension="fasta").output[0]
  job_ids = glob_wildcards(os.path.join(parallel_dir, "{index}.fasta")).index
  completed_files = expand("parallel/{task}/{index}.out",index=job_ids, task=wildcards["task"])
  return completed_files


rule annotation_merge_fasta:
  input:
    parallel_annotation_fasta
  output:
    "annotations_fasta/{task}.out"
  log:
    "logs/log_{task}_merge.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    echo input is {input}
    cat {input} > {output}
    """


rule annotation_merge_orfs:
  input:
    parallel_annotation_orfs
  output:
    "annotations_orfs/{task}.out"
  log:
    "logs/log_{task}_merge.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    echo input is {input}
    cat {input} > {output}
    """


rule pfam_parallel:
  input:
    "parallel/annotation_orfs/{index}.orfs",
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
    hmmscan --cpu {threads} --tblout {output} {input[1]} {input[0]}
    """


rule sprot_blastp_parallel:
  input:
    "parallel/annotation_orfs/{index}.orfs",
    "db/uniprot_sprot.fasta"
  output:
    "parallel/sprot_blastp/{index}.out"
  log:
    "logs/log_sprot_blastp_{index}.txt"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastp -query {input[0]} -db {input[1]} -num_threads {threads} -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output}
    """


rule sprot_blastx_parallel:
  input:
    "parallel/annotation_fasta/{index}.fasta",
    "db/uniprot_sprot.fasta"
  output:
    "parallel/sprot_blastx/{index}.out"
  log:
    "logs/log_sprot_blastx_{index}.txt"
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastx -query {input[0]} -db {input[1]} -num_threads {threads} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output}
    """


rule tmhmm_parallel:
  input:
    "parallel/annotation_orfs/{index}.orfs"
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
    tmhmm --short < {input} > {output}
    """


rule signalp_parallel:
  input:
    "parallel/annotation_orfs/{index}.orfs"
  output:
    "parallel/signalp/{index}.out"
  log:
    "logs/log_signalp_{index}.txt"
  params:
    memory="2"
  threads:
    1
  shell:
    """
    signalp -t {SIGNALP_ORGANISM} -f short {input} > {output}
    """


rule kallisto:
  input:
    samples="samples.txt",
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
    {TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts {input.transcriptome} {STRAND_SPECIFIC} --seqType fq --samples_file {input.samples} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input.gene_trans_map} 2>&1 > {log}
    {TRINITY_HOME}/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input.gene_trans_map} */abundance.tsv 2>&1 >> {log}
    if [ -f kallisto.isoform.TMM.EXPR.matrix ]; then
      cp kallisto.isoform.TMM.EXPR.matrix {output[0]}
    elif [ -f kallisto.isoform.TPM.not_cross_norm ]; then
      cp kallisto.isoform.TPM.not_cross_norm {output[0]} 
    else
      echo Neither kallisto.isoform.TMM.EXPR.matrix or kallisto.isoform.TPM.not_cross_norm were produced
      exit 1
    fi
    if [ -f kallisto.gene.TMM.EXPR.matrix ]; then
      cp kallisto.gene.TMM.EXPR.matrix {output[1]}
    elif [ -f kallisto.gene.TPM.not_cross_norm ]; then
      cp kallisto.gene.TPM.not_cross_norm {output[1]} 
    else
      echo Neither kallisto.gene.TMM.EXPR.matrix or kallisto.gene.TPM.not_cross_norm were produced
      exit 1
    fi
    """


rule annotated_fasta:
  input:
    "transcriptome.fasta",
    "transcriptome.pep",
    "transcriptome_expression_isoform.tsv",
    "annotations_fasta/sprot_blastx.out",
    "annotations_orfs/sprot_blastp.out",
    "annotations_orfs/pfam.out",
    "annotations_orfs/tmhmm.out",
    "annotations_orfs/signalp.out"
  output:
    "transcriptome_annotated.fasta",
    "transcriptome_annotated.pep"
  log:
    "logs/log_annotated_fasta.txt"
  params:
    memory="4"
  threads:
    1
  run:
    ## Annotation map: transcript id -> description
    transcript_annotations = {}
    protein_annotations = {}

    ## Load kallisto results
    print ("Loading expression values from", input[2])
    with open(input[2]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      columns = next(csv_reader)
      for row in csv_reader:
        annotation = "TPM:"
        for i in range(1, len(columns)):
          annotation += " " + columns[i] + "=" + str(row[i])
        transcript_annotations[row[0]] = transcript_annotations.get(row[0], "") + "<br>" + annotation

    ## Load blastx results
    print ("Loading blastx results from", input[3])
    with open(input[3]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      for row in csv_reader:
        if (len(row) < 13): continue
        annotation = "blastx: sp" + row[1] + " " + row[12] + "; e=" + str(row[10])
        transcript_annotations[row[0]] = transcript_annotations.get(row[0], "") + "<br>" + annotation

    ## Load blastp results
    print ("Loading blastp results from", input[4])
    with open(input[4]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      for row in csv_reader:
        if (len(row) < 13): continue
        annotation = "blastp: sp" + row[1] + " " + row[12] + "; e=" + str(row[10])
        protein_annotations[row[0]] = protein_annotations.get(row[0], "") + "<br>" + annotation

    ## Load pfam results
    print ("Loading pfam predictions from", input[5])
    with open(input[5]) as input_handle:
      for line in input_handle:
        if (line.startswith("#")): continue
        row = re.split(" +", line, 22)
        if (len(row) < 23): continue
        annotation = "pfam: " + row[1] + " " + row[22] + "; e=" + str(row[6])
        protein_annotations[row[3]] = protein_annotations.get(row[3], "") + "<br>" + annotation

    ## Load tmhmm results
    print ("Loading tmhmm predictions from", input[6])
    with open(input[6]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      for row in csv_reader:
        if (len(row) < 6): continue
        annotation = "tmhmm: " + row[2] + "; " + row[3] + "; " + row[4] + "; " + row[5]
        protein_annotations[row[0]] = protein_annotations.get(row[0], "") + "<br>" + annotation
    
    ## Load signalp results
    print ("Loading signalp predictions from", input[7])
    with open(input[7]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      for row in csv_reader:
        if (len(row) < 9): continue
        annotation = "signalp: " + str(row[5]) 
        protein_annotations[row[0]] = protein_annotations.get(row[0], "") + "<br>" + annotation
    
    ## Do the work
    print ("Annotating FASTA file", input[0], "to", output[0])
    with open(input[0], "r") as input_fasta_handle, open(output[0], "w") as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        record.description = transcript_annotations.get(record.id, "")
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")
    
    print ("Annotating FASTA file", input[1], "to", output[1])
    with open(input[1], "r") as input_fasta_handle, open(output[1], "w") as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = re.sub("\.p[0-9]+$", "", record.id)
        record.description = transcript_annotations.get(transcript_id, "") + protein_annotations.get(record.id, "")
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")


rule transcriptome_TPM_blast_table:
  input:
    "transcriptome_expression_isoform.tsv",
    "annotations_fasta/sprot_blastx.out"
  output:
    "transcriptome_TPM_blast.csv"
  log:
    "logs/log_transcriptome_TPM_blast_table.txt"
  params:
    memory="2"
  threads:
    1
  run:
    ## Annotation map: transcript id -> description
    blastx_annotations = {}

    ## Load blastx results
    print ("Loading blastx results from", input[1])
    with open(input[1]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter="\t")
      for row in csv_reader:
        if (len(row) < 13): continue
        blastx_annotations[row[0]] = row[12] + "; e=" + str(row[10])
    
    with open(input[0], "r") as input_csv_handle, open(output[0], "w") as output_csv_handle:
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
    wget --directory-prefix db "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz" 2>&1 > {log}
    gunzip db/uniprot_sprot.fasta.gz 2>&1 >> {log}
    makeblastdb -in db/uniprot_sprot.fasta -dbtype prot 2>&1 >> {log}
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
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz" 2>&1 > {log}
    gunzip db/Pfam-A.hmm.gz 2>&1 >> {log}
    hmmpress db/Pfam-A.hmm 2>&1 >> {log}
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
    wget --directory-prefix db "http://eggnogdb.embl.de/download/latest/data/NOG/NOG.annotations.tsv.gz" 2>&1 > {log}
    gunzip db/NOG.annotations.tsv.gz 2>&1 >> {log}
    """


