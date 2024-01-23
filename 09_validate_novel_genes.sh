#!/usr/bin/env bash

#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=500M
#SBATCH --time=02:00:00
#SBATCH --partition=pall
#SBATCH --job-name=bedtools
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_bedtools_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_bedtools_%j.e

GTF_DIR=/data/users/lfrei/rna_seq/assembly_results/all_merged.gtf
BED_DIR=/data/users/lfrei/rna_seq/bed_files
OUTPUT_DIR=/data/users/lfrei/rna_seq/integrative_analysis
REFERENCE_DIR=/data/courses/rnaseq_course/lncRNAs/Project1/references

mkdir -p $BED_DIR
mkdir -p $OUTPUT_DIR

#load modules
module load UHTS/Analysis/BEDTools/2.29.2
module load SequenceAnalysis/GenePrediction/cpat/1.2.4

#create bed6 file from the merged gtf (with chromosome, start, end, gene_id, score, strand )
#-1 because GTF files are 1-based, but bed files are 0-based
grep 'chr' $GTF_DIR | awk 'BEGIN {OFS="\t"} $3=="transcript" {print $1, $4-1, $5-1, $12, $6, $7}' > $BED_DIR/all.bed

#split all.bed into annotated and novel features
grep 'ENS' $BED_DIR/all.bed > $BED_DIR/annotated.bed
grep 'MSTR' $BED_DIR/all.bed > $BED_DIR/novel.bed

#create bed for 3' and 5' ends from novel features
awk 'BEGIN {OFS="\t"} $6=="+" {print $1, $2-50, $2+50, $4, $5, $6} $6=="-" {print $1, $3-50, $3+50, $4, $5, $6}' $BED_DIR/novel.bed > $BED_DIR/novel_5prime.bed
awk 'BEGIN {OFS="\t"} $6=="+" {print $1, $3-50, $3+50, $4, $5, $6} $6=="-" {print $1, $2-50, $2+50, $4, $5, $6}' $BED_DIR/novel.bed > $BED_DIR/novel_3prime.bed

#look for intergenic transcripts
bedtools intersect -v -wa -a $BED_DIR/novel.bed -b $BED_DIR/annotated.bed > $OUTPUT_DIR/intergenic.txt

#validate 5' annotations (transcription start site)
bedtools intersect -s -wa -a $BED_DIR/novel_5prime.bed -b $REFERENCE_DIR/refTSS_v4.1_human_coordinate.hg38.bed > $OUTPUT_DIR/annotation_5prime.txt

#validate 3' annotations (poly A)
bedtools intersect -s -wa -a $BED_DIR/novel_3prime.bed -b $REFERENCE_DIR/atlas.clusters.2.0.GRCh38.96.bed > $OUTPUT_DIR/annotation_3prime.txt

#asses protein-coding potential with cpat
# get fasta of bed file with novel transcripts
bedtools getfasta -s -name -fi $REFERENCE_DIR/GRCh38.genome.fa -bed $BED_DIR/novel.bed -fo $OUTPUT_DIR/novel_transcripts.fa

# run cpat
cpat.py -x $REFERENCE_DIR/Human_Hexamer.tsv -d $REFERENCE_DIR/Human_logitModel.RData -g $OUTPUT_DIR/novel_transcripts.fa -o $OUTPUT_DIR/protein_coding_potential