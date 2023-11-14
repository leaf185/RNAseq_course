#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pall
#SBATCH --job-name=assembly_stringtie
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,error
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_hisat2_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_hisat2_%j.e
#SBATCH --array=0-5

#load module
module load UHTS/Aligner/stringtie/1.3.3b

REF_ANNOTATION=/data/users/lfrei/rna_seq/reference/gencode_v44_chr_annotation_all.gff3.gz
OUTPUT_DIR=/data/users/lfrei/rna_seq/assembly_results
INPUT_DIR=/data/users/lfrei/rna_seq/alignment_results

mkdir -p $OUTPUT_DIR

#create file array
cd $INPUT_DIR
FILES=( $(ls | grep '.bam$') )-p
sample_name="${FILES[$SLURM_ARRAY_TASK_ID]%%.bam}"


#run stringtie on the 6 samples
stringtie -o ${sample_name}.gtf \
    -p 4 --rf\
    -G $REF_ANNOTATION_DIR/gencode_v44_chr_annotation_all.gff3\
    ${FILES[$SLURM_ARRAY_TASK_ID]}
