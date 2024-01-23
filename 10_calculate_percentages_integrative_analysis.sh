#!/usr/bin/env bash

#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=100M
#SBATCH --time=01:00:00
#SBATCH --partition=pall
#SBATCH --job-name=integrative_analysis_percentages
#SBATCH --mail-user=lea.frei@students.unibe.ch
#SBATCH --mail-type=begin,end,fail
#SBATCH --output=/data/users/lfrei/rna_seq/output/output_integrative_analysis_percentages_%j.o
#SBATCH --error=/data/users/lfrei/rna_seq/errors/error_integrative_analysis_percentages_%j.e

BED_DIR=/data/users/lfrei/rna_seq/bed_files
WORKING_DIR=/data/users/lfrei/rna_seq/integrative_analysis

novel_transcripts=$(wc -l < $BED_DIR/novel.bed)
correct_5prime=$(wc -l < $WORKING_DIR/annotation_5prime.txt)
percentage_5prime=$(echo "scale=2; ($correct_5prime/$novel_transcripts) * 100" | bc)
correct_3prime=$(wc -l < $WORKING_DIR/annotation_3prime.txt)
percentage_3prime=$(echo "scale=2; ($correct_3prime/$novel_transcripts) * 100" | bc)
intergenic=$(wc -l < $WORKING_DIR/intergenic.txt)
percentage_intergenic=$(echo "scale=2; ($intergenic/$novel_transcripts) * 100" | bc)

echo "Number of novel transcripts: $novel_transcripts" > $WORKING_DIR/integrative_analysis_percentages.txt
echo "Correct 5prime annotation: $correct_5prime, $percentage_5prime %" >> $WORKING_DIR/integrative_analysis_percentages.txt
echo "Correct 3prime annotation: $correct_3prime, $percentage_3prime %" >> $WORKING_DIR/integrative_analysis_percentages.txt
echo "Intergenic transcripts: $intergenic, $percentage_intergenic %" >> $WORKING_DIR/integrative_analysis_percentages.txt