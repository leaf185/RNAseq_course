#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pall
#SBATCH --job-name=merge_stringtie
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,error
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_stringtie_merge_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_stringtie_merge_%j.e

#load module
module load UHTS/Aligner/stringtie/1.3.3b

REF_ANNOTATION=/data/courses/rnaseq_course/lncRNAs/Project1/references/gencode.v44.annotation.gtf
INPUT_DIR=/data/users/lfrei/rna_seq/assembly_results
FILES=list_gtf_files.txt #in same folder as input .gtf files

cd $INPUT_DIR

#merge
stringtie --merge \
    -G $REF_ANNOTATION \
    -o all_merged.gtf \
    $FILES
