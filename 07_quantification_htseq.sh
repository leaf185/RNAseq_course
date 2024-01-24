#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G
#SBATCH --time=03:00:00
#SBATCH --partition=pall
#SBATCH --job-name=quantification_htseq
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_htseq_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_htseq_%j.e
#SBATCH --array=0-5

#load module
module load UHTS/Analysis/HTSeq/0.9.1

INPUT_DIR=/data/users/lfrei/rna_seq/alignment_results
OUTPUT_DIR=/data/users/lfrei/rna_seq/quantification_results
GTF_REFERENCE_DIR=/data/users/lfrei/rna_seq/assembly_results

mkdir -p $OUTPUT_DIR

#create file array
cd $INPUT_DIR
FILES=( $(ls | grep '.sam$') )
sample_name="${FILES[$SLURM_ARRAY_TASK_ID]%%.sam}"

#run htseq-count on each of the 6 samples to count by gene
htseq-count --mode=union \
    --format=sam \
    --nonunique=none \
    --stranded=reverse \
    --type=exon \
    --idattr=gene_id  \
    ${FILES[$SLURM_ARRAY_TASK_ID]} \
    $GTF_REFERENCE_DIR/all_merged.gtf > $OUTPUT_DIR/${sample_name}_by_gene.tsv

#run htseq-count on each of the 6 samples to count by transcript
htseq-count --mode=union \
    --format=sam \
    --nonunique=none \
    --stranded=reverse \
    --type=exon \
    --idattr=transcript_id \
    ${FILES[$SLURM_ARRAY_TASK_ID]} \
    $GTF_REFERENCE_DIR/all_merged.gtf > $OUTPUT_DIR/${sample_name}_by_transcript.tsv