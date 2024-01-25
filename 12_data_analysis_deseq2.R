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

setwd("C:/Users/freil/Documents/1_Studium/23_HS/RNA Sequencing/Project_lncRNA/tsv_data")

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

#add gene names etc. (code from DESeq2 guide)
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

####some stuff from cancer genomics----------------------------------
library(tibble)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(GSVA)
library(enrichplot)
library(ggpubr)
library(ggrepel)

#create matrix of normalized counts
normalized_counts <- counts(dds_genes, normalized=TRUE)
# - make normalised count matrix into dataframe
normalized_counts_df <- as.data.frame(counts(dds_genes, normalized=TRUE))

#make dataframes of only parental/only samples, add a column specifying condition, then merge them back together
normalized_counts_parental <- normalized_counts[,4:6]
normalized_counts_holoclonal <- normalized_counts[,1:3]
normalized_counts_parental <- cbind(normalized_counts_parental, condition = "parental")
colnames(normalized_counts_parental) <- c("1", "2", "3", "condition")
normalized_counts_holoclonal <- cbind(normalized_counts_holoclonal, condition = "holoclonal")
colnames(normalized_counts_holoclonal) <- c("1", "2", "3", "condition")
normalized_counts_all <- rbind(normalized_counts_parental, normalized_counts_holoclonal)
normalized_counts_all <- as.data.frame(normalized_counts_all)

normalized_counts_all[which(rownames(normalized_counts_all)=="MSTRG.26657"),] %>%
  ggplot(normalized_counts_all, aes(fill=condition)) +
    geom_boxplot()

#GSEA ---------------------------------------------

results_genes <- as.data.frame(results_genes)

## generate gene list
list_gene = results_genes$stat
names(list_gene) = as.character(results_genes$entrez)
list_gene = sort(list_gene, decreasing = TRUE)


## import geneset lists
hallmarks <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene) 

reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME") %>% 
  dplyr::select(gs_name, entrez_gene) 

go <- msigdbr(species = "Homo sapiens", category = "C5") %>% 
  dplyr::select(gs_name, entrez_gene) 

go_bp <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:BP") %>% 
  dplyr::select(gs_name, entrez_gene)

go_cc <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:CC") %>% 
  dplyr::select(gs_name, entrez_gene)

go_mf <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "GO:MF") %>% 
  dplyr::select(gs_name, entrez_gene)

kegg <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG") %>% 
  dplyr::select(gs_name, entrez_gene)



############################################## GSEA #################################################

# perform GSEA using GSEA function
em2 <- GSEA(list_gene, TERM2GENE = hallmarks, pvalueCutoff = 1, eps = 0)  ## pvalue cut off 1 is important 

## put results of GSEA into dataframe
em2_df <- as.data.frame(em2)
# categorise pathways by significance level
em2_df <- em2_df %>%
  mutate("significance"=ifelse(p.adjust<=0.05 & p.adjust>0.01, "Sig p.adj <0.05",
                               ifelse(p.adjust<=0.01, "Sig p.adj <0.01", "ns")))

## bar plot
plot <- em2_df %>%
  ggplot(aes(x=NES, y=reorder(Description, NES), fill=factor(significance))) +
  geom_bar(stat = "identity", alpha=0.75) +
  scale_fill_manual(values=c("blue",  "red", "pink")) +
  labs(title = "GSEA: Hallmarks - holoclonal vs parental",
       y="Hallmarks pathways",
       x="Normalised enrichment score (NES)") +
  theme_bw() +
  theme(axis.text = element_text(size = 8))   +
  guides(fill=guide_legend(title="Significance"))
plot

############################################## Edge plots #################################################

## get specific pathway ID
em2@result[["ID"]]
# generate plot with one pathway ID and display pvalues
gseaplot2(em2, geneSetID = 1, pvalue_table = TRUE)
# generate plot with >1 pathway ID and display pvalues
gseaplot2(em2, geneSetID = c(1:3), subplots = 1:2, pvalue_table = TRUE)

## visualise GSEA
em2_df <- as.data.frame(em2)
# categorise by significance level
em2_df <- em2_df %>%
  mutate("significance"=ifelse(p.adjust<=0.05 & p.adjust>0.01, "Sig p.adj <0.05",
                               ifelse(p.adjust<=0.01, "Sig p.adj <0.01", "ns")))

############################################## OVER REPRESENTATION ANALYSIS AND VISUALISATION #################################################

## filter signficantly enriched genes 
res_ora_up <- results_genes %>%
  filter(padj<=0.1 & log2FoldChange > 0.5)

res_ora_up <- as.data.frame(res_ora_up)

## generate gene list
list_gene = res_ora_up$stat
names(list_gene) = as.character(res_ora_up$entrez)
list_gene = sort(list_gene, decreasing = TRUE)

# this is specific to ORA analysis
gene <- names(list_gene)

#perform ORA using "enricher" function
em <- enricher(gene, TERM2GENE=hallmarks, pvalueCutoff = 0.1, 
               pAdjustMethod = "BH")

## visualise ORA using upsetplot
upsetplot(em)