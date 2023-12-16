#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=count_expressed_genes
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,error
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_expressed_genes_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_expressed_genes_%j.e

INPUT_DIR=/data/users/lfrei/rna_seq/quantification_results
OUTPUT_DIR=/data/users/lfrei/rna_seq

cd $INPUT_DIR

for file in *; do
    awk -v file=$file 'BEGIN {count=0} ($2>0) {count++} END {print file, "total", count}' $file >> $OUTPUT_DIR/expressed_genes_from_sam.txt
    grep 'ENS' $file | awk -v file=$file 'BEGIN {count=0} ($2>0) {count++} END {print file, "novel", count}' >> $OUTPUT_DIR/expressed_genes_from_sam.txt;
done