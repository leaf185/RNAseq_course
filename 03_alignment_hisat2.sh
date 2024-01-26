#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=50G
#SBATCH --time=04:00:00
#SBATCH --partition=pall
#SBATCH --job-name=alignment_hisat2
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_hisat2_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_hisat2_%j.e
#SBATCH --array=0-5

#load module
module load UHTS/Aligner/hisat/2.2.1
module load UHTS/Analysis/samtools/1.10

REFERENCE_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/references
OUTPUT_DIR=/data/users/lfrei/rna_seq/alignment_results
READS_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
HISAT2_INDEXES=/data/users/lfrei/rna_seq/reference_indexed

mkdir -p $OUTPUT_DIR

#index file of reference genome
#options: -f input is fastq, -p to use parallel threads
hisat2-build -f -p 2  $REFERENCE_DIR/GRCh38.genome.fa reference_indexed/GRCh38_index

#create file arrays for R1 and R2 reads
cd $READS_DIR
files_R1=( $(ls | grep '^[1P].*_R1_') )
files_R2=( $(ls | grep '^[1P].*_R2_') )
sample_name="${files_R1[$SLURM_ARRAY_TASK_ID]%%_R1*}"

#run alignment
cd $READS_DIR
hisat2 --rna-strandness RF -p 2 \
        -x $HISAT2_INDEXES/GRCh38_index \
        -1 "${files_R1[$SLURM_ARRAY_TASK_ID]}" \
        -2 "${files_R2[$SLURM_ARRAY_TASK_ID]}" \
        -S "$OUTPUT_DIR/$sample_name.sam"

#convert SAM to BAM and sort BAM file
samtools view -bS "$OUTPUT_DIR/$sample_name.sam" | samtools sort -o "$OUTPUT_DIR/$sample_name.bam"

#index sorted BAM file
samtools index "$OUTPUT_DIR/$sample_name.bam"
