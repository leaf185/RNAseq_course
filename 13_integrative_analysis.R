library("tidyverse")
library("ggplot2")
library("stringr")

#run the output script of CPAT to obtain the final table with protein coding probability
load("Human_logitModel.RData")
test <- read.table(file="protein_coding_potential.dat",sep="\t",col.names=c("ID","mRNA","ORF","Fickett","Hexamer"))
test$prob <- predict(mylogit,newdata=test,type="response")
attach(test)
output <- cbind("mRNA_size"=mRNA,"ORF_size"=ORF,"Fickett_score"=Fickett,"Hexamer_score"=Hexamer,"coding_prob"=test$prob)
write.table(output,file="protein_coding_potential.txt",quote=F,sep="\t",row.names=ID)

#load all the datasets
gene_table_data <- read.table("gene_table_filtered.txt", header = T, sep = ";")
de_results_data <- read.table("gene_expression_table.csv", header = T, sep = ",")
five_prime_data <- read.table("annotation_5prime.txt", sep = "\t")
three_prime_data <- read.table("annotation_3prime.txt", sep = "\t")
intergenic_data <- read.table("intergenic.txt", sep = "\t")
protein_coding_potential_data <- read.table("protein_coding_potential.txt", sep = "\t")

#calculate number of protein coding novel transcripts (at cutoff 0.364)
protein_coding_transcripts <- sum(protein_coding_potential_data$coding_prob >= 0.364)

#select relevant columns and remove ; after transcript ids
gene_table <- dplyr::select(gene_table_data, gene_id, transcript_id, biotype)

de_results <- de_results_data %>%
  dplyr::select(X, baseMean, log2FoldChange, padj, ensembl, entrez, hgnc_symbol) %>%
  rename("gene_id" = X)
    
protein_coding_potential <- protein_coding_potential_data %>%
  mutate(transcript_id = sapply(strsplit(rownames(protein_coding_potential_data), split="\\;"), "[", 1 )) %>%
  dplyr::select(transcript_id, coding_prob)

five_prime <- five_prime_data %>%
  mutate(transcript_id = sapply(strsplit(V4, split="\\;"), "[", 1)) %>%
  dplyr::select(transcript_id)

three_prime <- three_prime_data %>%
  mutate(transcript_id = sapply(strsplit(V4, split="\\;"), "[", 1)) %>%
  dplyr::select(transcript_id)

intergenic <- intergenic_data %>%
  mutate(transcript_id = sapply(strsplit(V4, split="\\;"), "[", 1)) %>%
  dplyr::select(transcript_id)


#add the additional information from DESeq to the gene_table¨
complete_table <- gene_table %>%
  left_join(de_results, by = "gene_id")

#filter for significantly different expression to find interesting genes in general
complete_table_sorted <- complete_table %>%
  filter(padj < 0.1) %>%
  arrange(baseMean) %>%
  arrange(desc(abs(log2FoldChange)))

#show top 20 candidates from this
head(complete_table_sorted, n = 20)

#subset for novel genes
novel_genes <- complete_table[grepl("MSTRG", complete_table$gene_id),] %>%
  dplyr::select(gene_id, transcript_id, baseMean, log2FoldChange, padj)

#add protein coding potential and columns with booleans for 3' and 5' annotation and intergenic
novel_genes_all_info <- novel_genes %>%
  mutate(transcript_id = str_trim(transcript_id)) %>%
  left_join(protein_coding_potential, by = "transcript_id") %>%
  mutate(correct_5prime = transcript_id %in% five_prime$transcript_id) %>%
  mutate(correct_3prime = transcript_id %in% three_prime$transcript_id) %>%
  mutate(intergenic = transcript_id %in% intergenic$transcript_id)

#sort novel genes to find promising new lncRNA therapeutic targets
novel_genes_sorted <- novel_genes_all_info %>%
  filter(padj < 0.1) %>%
  rowwise() %>%
  mutate(sum_for_rank = sum(correct_5prime, correct_3prime, intergenic)) %>%
  arrange(desc(sum_for_rank),
          desc(abs(log2FoldChange)), desc(baseMean), coding_prob)

write.csv(novel_genes_sorted, file = "novel_genes_ranked.csv")
