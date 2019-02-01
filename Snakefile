# Note: must avoid spaces in filenames in samples.txt

import os
import shutil
import re
import csv
import Bio.SeqIO
import Bio.Alphabet

shell.prefix("set -euo pipefail;")


# Note: there is another signalp in /usr/local/bin, but that one will crash with > 10000 sequences
SIGNALP_ORGANISM="euk"


TRINITY_HOME=os.path.dirname(shutil.which("Trinity"))

STRAND_SPECIFIC = "" # --SS_lib_type RF

TRINITY_PARAMS = " --seqType fq"
TRINITY_PARAMS += " --trimmomatic --quality_trimming_params \"ILLUMINACLIP:" + TRINITY_HOME + "/trinity-plugins/Trimmomatic/adapters/TruSeq3-PE-2.fa:3:30:10 SLIDINGWINDOW:4:20 LEADING:20 TRAILING:20 MINLEN:25\""
TRINITY_PARAMS += " " + STRAND_SPECIFIC
TRINITY_PARAMS += " --min_glue 2"
TRINITY_PARAMS += " --min_kmer_cov 3"
#TRINITY_PARAMS += " --full_cleanup"
#TRINITY_PARAMS += " --no_normalize_reads"

# https://github.com/griffithlab/rnaseq_tutorial/wiki/Trinity-Assembly-And-Analysis


rule all:
  input:
    'samples.txt',
    'transcriptome_stats.txt',
    'transcriptome_exN50.tsv.plot.pdf',
    'transcriptome_annotated.fasta',
    'proteome_annotated.fasta',
    'transcriptome_TPM_blast.csv'


rule clean:
  shell:
    """
    rm -rf trinity_* tmp* log* TMHMM* kallisto* transcriptome* proteome* transrate_results* pfam* sprot* signalp* tmhmm* db*
    if [ -f samples.txt ]; then
      cut -f 2 < samples.txt | xargs --no-run-if-empty rm -rf
    fi
    """


rule trinity_inchworm_chrysalis:
  input:
    'samples.txt'
  output:
    'trinity_out_dir/recursive_trinity.cmds'
  log:
    'logs/log_trinity_inchworm_chrysalis.txt'
  params:
    memory="200"
  threads:
    16
  shell:
    """
    Trinity --no_distributed_trinity_exec --max_memory {params.memory}G --CPU {threads} --samples_file {input} {TRINITY_PARAMS}
    """


rule trinity_butterfly_split:
  input:
    'trinity_out_dir/recursive_trinity.cmds'
  output:
    dynamic('trinity_out_dir/parallel/job{job_index}')
  #log:
  #    'logs/log_trinity_stage2.txt'
  params:
    memory="10"
  threads:
    1
  shell:
    """
    split -l 100 {input} trinity_out_dir/parallel/job
    """


rule trinity_butterfly_parallel:
  input:
    'trinity_out_dir/parallel/job{job_index}'
  output:
    'trinity_out_dir/parallel/completed{job_index}'
  #log:
  #    'logs/log_trinity_stage3.txt'
  params:
    memory="10"
  threads:
    1
  shell:
    """
    bash {input}
    cat {input} > {output}
    """


rule trinity_butterfly_merge:
  input:
    'samples.txt',
    'trinity_out_dir/recursive_trinity.cmds',
    dynamic('trinity_out_dir/parallel/completed{job_index}'),
  output:
    'trinity_out_dir/recursive_trinity.cmds.completed',
    'transcriptome.fasta',
    'transcriptome.gene_trans_map'
  log:
    'logs/log_trinity_butterfly_merge.txt'
  params:
    memory="200"
  threads:
    16
  shell:
    """
    cp -n {input[1]} {output[0]}
    Trinity --max_memory {params.memory}G --CPU {threads} --samples_file {input[0]} {TRINITY_PARAMS}
    mv trinity_out_dir/Trinity.fasta {output[1]}
    mv trinity_out_dir/Trinity.fasta.gene_trans_map {output[2]}
    """


rule trinity_stats:
  input:
    'transcriptome.fasta',
    'transcriptome_expression_isoform.tsv'
  output:
    'transcriptome_stats.txt',
    'transcriptome_exN50.tsv',
    'transcriptome_exN50.tsv.plot.pdf'
  log:
    'logs/log_trinity_exN50.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    {TRINITY_HOME}/util/TrinityStats.pl {input[0]} > {output[0]}
    {TRINITY_HOME}/util/misc/contig_ExN50_statistic.pl {input[1]} {input[0]} > {output[1]}
    {TRINITY_HOME}/util/misc/plot_ExN50_statistic.Rscript {output[1]}
    """


rule transrate:
  input:
    'samples.txt',
    'transcriptome.fasta'
  output:
    'transrate_results/transcriptome/contigs.csv'
  log:
    'logs/log_transrate.txt'
  params:
    memory="2"
  threads:
    12
  shell:
    """
    LEFT=`cut -f 3 < {input[0]} | tr '\n' ',' | sed 's/,*$//g'`
    RIGHT=`cut -f 4 < {input[0]} | tr '\n' ',' | sed 's/,*$//g'`
    transrate --threads {threads} --assembly={input[1]} --left=$LEFT --right=$RIGHT
    """


rule transdecoder:
  input:
    'transcriptome.fasta'
  output:
    'proteome.fasta'
  log:
    'logs/log_transdecoder.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    rm -rf {input}.transdecoder_dir
    TransDecoder.LongOrfs -t {input}
    mv {input}.transdecoder_dir/longest_orfs.pep {output}
    """


rule align_reads:
  input:
    'samples.txt',
    'transcriptome.fasta',
    'transcriptome.gene_trans_map'
  output:
    'RSEM_out.gene.TMM.EXPR.matrix'
  log:
    'logs/log_bowtie2.txt'
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
    'samples.txt',
    'kallisto.gene.counts.matrix'
  output:
    'edgeR_trans'
  log:
    'logs/log_trinity_DE.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    {TRINITY_HOME}/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix {input[1]} --samples_file {input[1]} --method edgeR --output {output}
    """


rule pfam_split:
  input:
    'proteome.fasta'
  output:
    dynamic('pfam/proteome{pfam_index}.fasta')
  #log:
  #    'logs/log_pep_split.txt'
  params:
    memory="2"
  threads:
    1
  run:
    with open(input[0], "r") as input_handle:
      parser = Bio.SeqIO.parse(input_handle, "fasta")
      count = 0
      for entry in parser:
        if (count % 1000 == 0):
          filename = output[0].replace("__snakemake_dynamic__", str(int(count / 1000)))
          handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        entry.seq = entry.seq.strip("*")
        Bio.SeqIO.write(entry, handle, "fasta")
        count += 1


rule pfam_parallel:
  input:
    'pfam/proteome{pfam_index}.fasta',
    'db/Pfam-A.hmm'
  output:
    'pfam/proteome{pfam_index}_pfam.out'
  log:
    'logs/log_pfam{pfam_index}.txt'
  params:
    memory="2"
  threads:
    2
  shell:
    """
    hmmscan --cpu {threads} --tblout {output} {input[1]} {input[0]}
    """


rule pfam_merge:
  input:
    dynamic('pfam/proteome{pfam_index}_pfam.out')
  output:
    'pfam/pfam.out'
  log:
    'logs/log_pfam_merge.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output}
    """


rule sprot_blastp_split:
  input:
    'proteome.fasta'
  output:
    dynamic('sprot_blastp/proteome{blastp_index}.fasta')
  #log:
  #    'logs/log_pep_split.txt'
  params:
    memory="2"
  threads:
    1
  run:
    with open(input[0], "r") as input_handle:
      parser = Bio.SeqIO.parse(input_handle, "fasta")
      count = 0
      for entry in parser:
        if (count % 1000 == 0):
          filename = output[0].replace("__snakemake_dynamic__", str(int(count / 1000)))
          handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        entry.seq = entry.seq.strip("*")
        Bio.SeqIO.write(entry, handle, "fasta")
        count += 1


rule sprot_blastp_parallel:
  input:
    'sprot_blastp/proteome{blastp_index}.fasta',
    'db/uniprot_sprot.fasta'
  output:
    'sprot_blastp/proteome{blastp_index}_sprot_blastp.out'
  log:
    'logs/log_sprot_blastp{blastp_index}.txt'
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastp -query {input[0]} -db {input[1]} -num_threads {threads} -evalue 1e-10 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output}
    """


rule sprot_blastp_merge:
  input:
    dynamic('sprot_blastp/proteome{blastp_index}_sprot_blastp.out')
  output:
    'sprot_blastp/sprot_blastp.out'
  log:
    'logs/log_sprot_blastp_merge.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output}
    """


rule sprot_blastx_split:
  input:
    'transcriptome.fasta'
  output:
    dynamic('sprot_blastx/proteome{blastx_index}.fasta')
  #log:
  #    'logs/log_fasta_split.txt'
  params:
    memory="2"
  threads:
    1
  run:
    with open(input[0], "r") as input_handle:
      parser = Bio.SeqIO.parse(input_handle, "fasta")
      count = 0
      for entry in parser:
        if (count % 1000 == 0):
          filename = output[0].replace("__snakemake_dynamic__", str(int(count / 1000)))
          handle = open(filename, "w");
        Bio.SeqIO.write(entry, handle, "fasta")
        count += 1


rule sprot_blastx_parallel:
  input:
    'sprot_blastx/proteome{blastx_index}.fasta',
    'db/uniprot_sprot.fasta'
  output:
    'sprot_blastx/proteome{blastx_index}_sprot_blastx.out'
  log:
    'logs/log_sprot_blastx{blastx_index}.txt'
  params:
    memory="4"
  threads:
    2
  shell:
    """
    blastx -query {input[0]} -db {input[1]} -num_threads {threads} -evalue 1e-6 -max_hsps 1 -max_target_seqs 1 -outfmt "6 std stitle" -out {output}
    """


rule sprot_blastx_merge:
  input:
    dynamic('sprot_blastx/proteome{blastx_index}_sprot_blastx.out')
  output:
    'sprot_blastx/sprot_blastx.out'
  log:
    'logs/log_sprot_blastx_merge.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output}
    """


rule tmhmm_split:
  input:
    'proteome.fasta'
  output:
    dynamic('tmhmm/proteome{tmhmm_index}.fasta')
  #log:
  #    'logs/log_pep_split.txt'
  params:
    memory="2"
  threads:
    1
  run:
    with open(input[0], "r") as input_handle:
      parser = Bio.SeqIO.parse(input_handle, "fasta")
      count = 0
      for entry in parser:
        if (count % 1000 == 0):
          filename = output[0].replace("__snakemake_dynamic__", str(int(count / 1000)))
          handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        entry.seq = entry.seq.strip("*")
        Bio.SeqIO.write(entry, handle, "fasta")
        count += 1


rule tmhmm_parallel:
  input:
    'tmhmm/proteome{tmhmm_index}.fasta',
  output:
    'tmhmm/proteome{tmhmm_index}_tmhmm.out'
  log:
    'logs/log_tmhmm{tmhmm_index}.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    tmhmm --short < {input} > {output}
    """


rule tmhmm_merge:
  input:
    dynamic('tmhmm/proteome{tmhmm_index}_tmhmm.out')
  output:
    'tmhmm/tmhmm.out'
  log:
    'logs/log_tmhmm_merge.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output}
    """


rule signalp_split:
  input:
    'proteome.fasta'
  output:
    dynamic('signalp/proteome{signalp_index}.fasta')
  #log:
  #    'logs/log_pep_split.txt'
  params:
    memory="2"
  threads:
    1
  run:
    with open(input[0], "r") as input_handle:
      parser = Bio.SeqIO.parse(input_handle, "fasta")
      count = 0
      for entry in parser:
        if (count % 1000 == 0):
          filename = output[0].replace("__snakemake_dynamic__", str(int(count / 1000)))
          handle = open(filename, "w");
        # Remove predicted stop codons, because some annotation tools do not like them (e.g. InterProScan) 
        entry.seq = entry.seq.strip("*")
        Bio.SeqIO.write(entry, handle, "fasta")
        count += 1


rule signalp_parallel:
  input:
    'signalp/proteome{signalp_index}.fasta',
  output:
    'signalp/proteome{signalp_index}_signalp.out'
  log:
    'logs/log_signalp{signalp_index}.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    signalp -t {SIGNALP_ORGANISM} -f short {input} > {output}
    """


rule signalp_merge:
  input:
    dynamic('signalp/proteome{signalp_index}_signalp.out')
  output:
    'signalp/signalp.out'
  log:
    'logs/log_signalp_merge.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    cat {input} > {output}
    """


rule kallisto:
  input:
    'samples.txt',
    'transcriptome.fasta',
    'transcriptome.gene_trans_map'
  output:
    'transcriptome_expression_isoform.tsv',
    'transcriptome_expression_gene.tsv'
  log:
    'logs/log_kallisto.txt'
  params:
    memory="2"
  threads:
    8
  shell:
    """
    {TRINITY_HOME}/util/align_and_estimate_abundance.pl --transcripts {input[1]} {STRAND_SPECIFIC} --seqType fq --samples_file {input[0]} --prep_reference --thread_count {threads} --est_method kallisto --gene_trans_map {input[2]}
    {TRINITY_HOME}/util/abundance_estimates_to_matrix.pl --est_method kallisto --name_sample_by_basedir --gene_trans_map {input[2]} */abundance.tsv
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
    'transcriptome.fasta',
    'proteome.fasta',
    'transcriptome_expression_isoform.tsv',
    'transrate_results/transcriptome/contigs.csv',
    'sprot_blastx/sprot_blastx.out',
    'sprot_blastp/sprot_blastp.out',
    'pfam/pfam.out',
    'tmhmm/tmhmm.out',
    'signalp/signalp.out'
  output:
    'transcriptome_annotated.fasta',
    'proteome_annotated.fasta'
  log:
    'logs/log_annotated_fasta.txt'
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
      csv_reader = csv.reader(input_handle, delimiter='\t')
      columns = next(csv_reader)
      for row in csv_reader:
        annotation = "TPM:"
        for i in range(1, len(columns)):
          annotation += " " + columns[i] + "=" + str(row[i])
        transcript_annotations[row[0]] = transcript_annotations.get(row[0], "") + "<br>" + annotation

    ## Load transrate results
    print ("Loading transrate results from", input[3])
    with open(input[3]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter=',')
      columns = next(csv_reader)
      for row in csv_reader:
        if (len(row) < 18): continue
        annotation = "transrate: " + columns[5] + "=" + str(row[5]) + " " + columns[7] + "=" + str(row[7])
        transcript_annotations[row[0]] = transcript_annotations.get(row[0], "") + "<br>" + annotation

    ## Load blastx results
    print ("Loading blastx results from", input[4])
    with open(input[4]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 13): continue
        annotation = "blastx: sp" + row[1] + " " + row[12] + "; e=" + str(row[10])
        transcript_annotations[row[0]] = transcript_annotations.get(row[0], "") + "<br>" + annotation

    ## Load blastp results
    print ("Loading blastp results from", input[5])
    with open(input[5]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 13): continue
        annotation = "blastp: sp" + row[1] + " " + row[12] + "; e=" + str(row[10])
        protein_annotations[row[0]] = protein_annotations.get(row[0], "") + "<br>" + annotation

    ## Load pfam results
    print ("Loading pfam predictions from", input[6])
    with open(input[6]) as input_handle:
      for line in input_handle:
        if (line.startswith("#")): continue
        row = re.split(" +", line, 22)
        if (len(row) < 23): continue
        annotation = "pfam: " + row[1] + " " + row[22] + "; e=" + str(row[6])
        protein_annotations[row[3]] = protein_annotations.get(row[3], "") + "<br>" + annotation

    ## Load tmhmm results
    print ("Loading tmhmm predictions from", input[7])
    with open(input[7]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 6): continue
        annotation = "tmhmm: " + row[2] + "; " + row[3] + "; " + row[4] + "; " + row[5]
        protein_annotations[row[0]] = protein_annotations.get(row[0], "") + "<br>" + annotation
    
    ## Load signalp results
    print ("Loading signalp predictions from", input[8])
    with open(input[8]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 9): continue
        annotation = "signalp: " + str(row[5]) 
        protein_annotations[row[0]] = protein_annotations.get(row[0], "") + "<br>" + annotation
    
    ## Do the work
    print ("Annotating FASTA file", input[0], "to", output[0])
    with open(input[0], 'r') as input_fasta_handle, open(output[0], 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        record.description = transcript_annotations.get(record.id, "")
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")
    
    print ("Annotating FASTA file", input[1], "to", output[1])
    with open(input[1], 'r') as input_fasta_handle, open(output[1], 'w') as output_fasta_handle:
      for record in Bio.SeqIO.parse(input_fasta_handle, "fasta"):
        transcript_id = re.sub("\.p[0-9]+$", "", record.id)
        record.description = transcript_annotations.get(transcript_id, "") + protein_annotations.get(record.id, "")
        Bio.SeqIO.write(record, output_fasta_handle, "fasta")


rule transcriptome_TPM_blast_table:
  input:
    'transcriptome_expression_isoform.tsv',
    'transrate_results/transcriptome/contigs.csv',
    'sprot_blastx/sprot_blastx.out'
  output:
    'transcriptome_TPM_blast.csv'
  log:
    'logs/log_transcriptome_TPM_blast_table.txt'
  params:
    memory="2"
  threads:
    1
  run:
    ## Annotation map: transcript id -> description
    transrate_pgood_annotations = {}
    blastx_annotations = {}

    ## Load transrate results
    print ("Loading transrate results from", input[1])
    with open(input[1]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter=',')
      transrate_columns = next(csv_reader)
      for row in csv_reader:
        if (len(row) < 18): continue
        transrate_pgood_annotations[row[0]] = row[5]

    ## Load blastx results
    print ("Loading blastx results from", input[2])
    with open(input[2]) as input_handle:
      csv_reader = csv.reader(input_handle, delimiter='\t')
      for row in csv_reader:
        if (len(row) < 13): continue
        blastx_annotations[row[0]] = row[12] + "; e=" + str(row[10])
    
    with open(input[0], 'r') as input_csv_handle, open(output[0], 'w') as output_csv_handle:
      csv_reader = csv.reader(input_csv_handle, delimiter='\t')
      csv_writer = csv.writer(output_csv_handle, delimiter=',')
      csv_columns = next(csv_reader)
      csv_columns[0] = "transcript"
      for i in range(1, len(csv_columns)):
        csv_columns[i] = "TPM(" + csv_columns[i] + ")"
      csv_columns.append("blastx")
      csv_columns.append("transrate pgood")
      csv_writer.writerow(csv_columns)
      for row in csv_reader:
        row.append(blastx_annotations.get(row[0], ""))
        row.append(transrate_pgood_annotations.get(row[0], ""))
        csv_writer.writerow(row)

 
rule download_sprot:
  output:
    'db/uniprot_sprot.fasta'
  log:
    'logs/log_download_sprot.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz"
    gunzip db/uniprot_sprot.fasta.gz
    makeblastdb -in db/uniprot_sprot.fasta -dbtype prot
    """


rule download_pfam:
  output:
    'db/Pfam-A.hmm'
  log:
    'logs/log_download_pfam.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/current_release/Pfam-A.hmm.gz"
    gunzip db/Pfam-A.hmm.gz
    hmmpress db/Pfam-A.hmm
    """


rule download_eggnog:
  output:
    'db/NOG.annotations.tsv'
  log:
    'logs/log_download_eggnog.txt'
  params:
    memory="2"
  threads:
    1
  shell:
    """
    wget --directory-prefix db "http://eggnogdb.embl.de/download/latest/data/NOG/NOG.annotations.tsv.gz"
    gunzip db/NOG.annotations.tsv.gz
    """


