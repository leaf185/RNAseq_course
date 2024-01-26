#install the package
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
BiocManager::install("EnhancedVolcano")
BiocManager::install("biomaRt")
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
library(biomaRt)
library(pheatmap)

#create sample_table to load htseq data 
sample_table_genes <- data.frame(sample = c("1_1", "1_2", "1_5", "P1", "P2", "P3"), 
                          filename = c("1_1_L3_by_gene.tsv", "1_2_L3_by_gene.tsv","1_5_L3_by_gene.tsv", "P1_L3_by_gene.tsv", "P2_L3_by_gene.tsv","P3_L3_by_gene.tsv"),
                          group = c("Holoclonal", "Holoclonal", "Holoclonal", "Parental", "Parental", "Parental"))
sample_table_transcripts <- data.frame(sample = c("1_1", "1_2", "1_5", "P1", "P2", "P3"), 
                                 filename = c("1_1_L3_by_transcript.tsv", "1_2_L3_by_transcript.tsv","1_5_L3_by_transcript.tsv", "P1_L3_by_transcript.tsv", "P2_L3_by_transcript.tsv","P3_L3_by_transcript.tsv"),
                                 group = c("Holoclonal", "Holoclonal", "Holoclonal", "Parental", "Parental", "Parental"))

#create DESeqDataSet from the HTSeq output
dds_genes <- DESeqDataSetFromHTSeqCount(sample_table_genes, design = ~ group)
dds_transcripts <- DESeqDataSetFromHTSeqCount(sample_table_transcripts, design = ~ group)

#make sure parental is the control level
dds_genes$group <- relevel(dds_genes$group, "Parental")
dds_transcripts$group <- relevel(dds_transcripts$group, "Parental")

#get the FPM (fragments per kilobase million)
dds_genes_fpm <- fpm(dds_genes)
dds_transcripts_fpm <- fpm(dds_transcripts)

#verify the read counts (fpm) add up to the expected number (1 million)
fpm_sums <- data.frame(sample = c("1_1", "1_2", "1_5", "P1", "P2", "P3"),
                       genes = colSums(dds_genes_fpm),
                       transcripts = colSums(dds_transcripts_fpm))

#differential expression analysis
dds_genes <- DESeq(dds_genes)
dds_transcripts <- DESeq(dds_transcripts)
results_genes <- results(dds_genes)
results_transcripts <- results(dds_transcripts)

#add gene names, biotypes, entrez gene id (code from DESeq2 guide)
results_genes$ensembl <- sapply( strsplit( rownames(results_genes), split="\\." ), "[", 1 )
results_transcripts$ensembl <- sapply( strsplit( rownames(results_transcripts), split="\\." ), "[", 1 )
ensembl = useMart( "ensembl", dataset = "hsapiens_gene_ensembl" )
genemap <- getBM( attributes = c("ensembl_gene_id", "entrezgene_id", "hgnc_symbol", "gene_biotype"),
                  filters = "ensembl_gene_id",
                  values = results_genes$ensembl,
                  mart = ensembl)
idx <- match( results_genes$ensembl, genemap$ensembl_gene_id )
results_genes$entrez <- genemap$entrezgene_id[ idx ]
results_genes$hgnc_symbol <- genemap$hgnc_symbol[ idx ]
results_genes$biotype <- genemap$gene_biotype [idx]
transcriptmap <- getBM(attributes = c("ensembl_transcript_id","entrezgene_id", "hgnc_symbol", "transcript_biotype"),
                       filters = "ensembl_transcript_id",
                       values = results_transcripts$ensembl,
                       mart = ensembl)
idy <- match (results_transcripts$ensembl, transcriptmap$ensembl_transcript_id)
results_transcripts$hgnc_symbol <- transcriptmap$hgnc_symbol[idy]
results_transcripts$biotype <- transcriptmap$transcript_biotype[idy]

#get number of significantly up- and downregulated genes
genes_significant_up <- as.data.frame(results_genes) %>%
  filter(padj < 0.1 & log2FoldChange > 0) %>%
  nrow()
genes_significant_down <- as.data.frame(results_genes) %>%
  filter(padj < 0.1 & log2FoldChange < 0) %>%
  nrow()
transcripts_significant_up <- as.data.frame(results_transcripts) %>%
  filter(padj < 0.1 & log2FoldChange > 0) %>%
  nrow()
transcripts_significant_down <- as.data.frame(results_transcripts) %>%
  filter(padj < 0.1 & log2FoldChange < 0) %>%
  nrow()

#save the expression tables
write.csv( as.data.frame(results_genes), file="gene_expression_table.csv" )
write.csv( as.data.frame(results_transcripts), file="transcript_expression_table.csv" )

#volcano plot for all genes/ all transcripts
EnhancedVolcano(results_genes,results_genes$hgnc_symbol, "log2FoldChange", "padj",
                xlim = c(-10,10)
                )
EnhancedVolcano(results_transcripts, results_transcripts$hgnc_symbol, "log2FoldChange", "padj",
                xlim = c(-10,10))

#volcano plot for lncRNAs (from gene expression)
EnhancedVolcano(results_genes[which(results_genes$biotype == "lncRNA"),], 
                results_genes$hgnc_symbol[which(results_genes$biotype == "lncRNA")], "log2FoldChange", "padj",
                xlim = c(-10,10))

#volcano plot for proteincoding genes (from gene expression)
EnhancedVolcano(results_genes[which(results_genes$biotype == "protein_coding"),],
                results_genes$hgnc_symbol[which(results_genes$biotype == "protein_coding")], "log2FoldChange", "padj",
                xlim = c(-10,10))

#volcano plot of all novel genes (from gene expression)
EnhancedVolcano(results_genes[grep("MSTRG", rownames(results_genes)), ],
                rownames(results_genes[grep("MSTRG", rownames(results_genes)), ]),
                "log2FoldChange", "padj", xlim = c(-10, 10))

#create heatmap from the fpm table for the top differential expressed genes
results_genes_significant <- results_genes[which(results_genes$padj < 0.1),]
results_genes_significant_ordered <- results_genes_significant[order(abs(results_genes_significant$log2FoldChange), decreasing = T),]
dds_genes_fpm_matrix <- as.matrix(dds_genes_fpm[which(rownames(dds_genes_fpm) %in% head(rownames(results_genes_significant_ordered))),])
pheatmap(dds_genes_fpm_matrix)
