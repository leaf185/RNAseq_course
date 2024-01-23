#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=100M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=count_genes
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_count_genes_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_count_genes_%j.e

INPUT_DIR=/data/users/lfrei/rna_seq/assembly_results/all_merged.gtf
OUTPUT_DIR=/data/users/lfrei/rna_seq

cd $INPUT_DIR

#number of exons
awk '$3=="exon" {exon_count++} END {print "exons:", exon_count}' $INPUT_DIR >> $OUTPUT_DIR/gene_count.txt

#number of transcripts
transcript_count=$(awk '$11=="transcript_id" {print $12}' $INPUT_DIR | sort | uniq | wc -l)
echo transcripts: $transcript_count >> $OUTPUT_DIR/gene_count.txt

#number of genes
gene_count=$(awk '$9=="gene_id" {print $10}' $INPUT_DIR | sort | uniq | wc -l)
echo genes: $gene_count >> $OUTPUT_DIR/gene_count.txt

#number of novel transcripts
novel_transcripts=$(awk '$11=="transcript_id" {print $12}' $INPUT_DIR | grep -v 'ENS' | sort | uniq | wc -l)
echo novel transcripts: $novel_transcripts >> $OUTPUT_DIR/gene_count.txt

#number of novel genes
novel_genes=$(awk '$9=="gene_id" {print $10}' $INPUT_DIR | grep -v 'ENS' | sort | uniq | wc -l)
echo novel genes: $novel_genes >> $OUTPUT_DIR/gene_count.txt 

#number of single exon genes 
awk '$3=="exon" {print $10}' $INPUT_DIR | sort | uniq -c | \
awk '$1==1 {single_exon_count++} END {print "single exon genes:", single_exon_count}' >> $OUTPUT_DIR/gene_count.txt

#number of single exon transcripts
awk '$3=="exon" {print $12}' $INPUT_DIR | sort | uniq -c | \
awk '$1==1 {single_exon_count++} END {print "single exon transcripts:", single_exon_count}' >> $OUTPUT_DIR/gene_count.txt