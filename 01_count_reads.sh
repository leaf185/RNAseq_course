#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=100M
#SBATCH --time=0:15:00
#SBATCH --partition=pall
#SBATCH --job-name=count_reads
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/reads.txt
#SBATCH --array=0-11

FASTQ_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
OUTPUT_DIR=/data/users/lfrei/rna_seq

cd $FASTQ_DIR
FILES=( $(ls | grep '^[1P]') ) #use all the files in the folder starting with either 1 or P

current_file="${FILES[$SLURM_ARRAY_TASK_ID]}"

count=$(zcat ${current_file} | awk 'BEGIN {count=0} /^@/ {++count} END {print count}')
echo "${current_file}, ${count}" >> $OUTPUT_DIR/reads.txt
