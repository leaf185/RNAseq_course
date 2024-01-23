#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=create_gene_table
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_create_table_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_create_table_%j.e

INPUT_DIR=/data/users/lfrei/rna_seq/assembly_results/all_merged.gtf
GTF_REFERENCE_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/references/gencode.v44.annotation.gtf
OUTPUT_DIR=/data/users/lfrei/rna_seq


#create table from the GENCODE GTF
echo -e "gene_id; transcript_id; biotype;" > $OUTPUT_DIR/gene_table_filtered.txt
grep -E 'lncRNA|protein_coding' $GTF_REFERENCE_DIR |
awk '$3=="transcript" {print $10, $12, $14}' >> $OUTPUT_DIR/gene_table_filtered.txt

#add novel genes from the merged GTF
#first grep step to keep only transcripts annotated on chromosomes (more reliable)
grep 'chr' $INPUT_DIR | awk '$3=="transcript" {print $10, $12, "novel;"}' | grep -v 'ENS' >> $OUTPUT_DIR/gene_table_filtered.txt


#create the table without removing any genes:
echo -e "gene_id; transcript_id; biotype;" > $OUTPUT_DIR/gene_table.txt
awk '$3=="transcript" {print $10, $12, $14}' $GTF_REFERENCE_DIR >> $OUTPUT_DIR/gene_table.txt

awk '$3=="transcript" {print $10, $12, "novel"}' $INPUT_DIR | grep -v 'ENS' >> $OUTPUT_DIR/gene_table.txt
