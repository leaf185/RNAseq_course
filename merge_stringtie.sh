#!/usr/bin/env bash

#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=8G
#SBATCH --time=04:00:00
#SBATCH --partition=pall
#SBATCH --job-name=merge_stringtie
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,error
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_hisat2_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/error/error_hisat2_%j.e

#load module
module load UHTS/Aligner/stringtie/1.3.3b

#merge
stringtie --merge list of gtf files -o outputfile.gtf\
    -G ref_annot-file