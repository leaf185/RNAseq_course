#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=100M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=quality_check
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,error
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_fastqc_%j.o
#SBATCH --array=0-11

module load UHTS/Quality_control/fastqc/0.11.9

READS_DIR=/data/courses/rnaseq_course/lncRNAs/fastq
OUTPUT_DIR=/data/users/lfrei/rna_seq/fastqc

mkdir -p $OUTPUT_DIR

cd $READS_DIR

FILES=( $(ls | grep '^[1P]') )

fastqc --outdir $OUTPUT_DIR ${FILES[$SLURM_ARRAY_TASK_ID]}