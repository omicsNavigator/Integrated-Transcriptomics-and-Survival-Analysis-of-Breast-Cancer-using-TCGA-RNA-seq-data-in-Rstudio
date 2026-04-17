#Data acquisition
library(TCGAbiolinks)
library(DESeq2)
query <- GDCquery(project = "TCGA-BRCA",data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",workflow.type = "STAR - Counts")

GDCdownload(query)
data <- GDCprepare(query)
data <- readRDS("tcga_brca.rds")
#saveRDS(data, "tcga_brca.rds") #saving the downloaded file to avoid redownload

#Extract expression matrix
counts <- assay(data)

#Define groups
colnames(colData(data)) #just seeing inside the coldata
metadata <- colData(data)

group <- ifelse(metadata$sample_type == "Primary Tumor","Tumor","Normal")
group <- factor(group)

table(group)
dds <- DESeqDataSetFromMatrix(countData = counts,colData = data.frame(condition = group),
                              design = ~ condition)

#Differential expression
dds <- DESeq(dds)
res <- results(dds)
summary(res)
res <- res[order(res$padj), ]
head(res)

#Significant genes
sigGenes <- res[!is.na(res$padj) & res$padj < 0.05, ]
nrow(sigGenes)

#Clean data by removing missing values
res_df <- as.data.frame(res)
res_df <- res_df[!is.na(res_df$pvalue), ]

#Volcano plot
library(ggplot2)
res_df$significant <- ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1,
                             "Significant","Not Significant")

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue),color = significant)) +
      geom_point(alpha = 0.5) +theme_minimal() +
      labs(title = "TCGA Breast Cancer Volcano Plot", x= "Log2 Fold Change",y="-log10(p-value)")

#Pathway Enrichment(GO + KEGG)
#install required packages
#BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))
#BiocManager::install("org.Hs.eg.db", lib = .libPaths()[1], force = TRUE)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)

#Prepare gene list & Convert ENSEMBL TO ENTREZ IDs
gene_list <- rownames(sigGenes)
head(rownames(sigGenes))
gene_list_clean <- gsub("\\..*", "", gene_list) 
head(gene_list_clean)

gene_df <- bitr(gene_list_clean,fromType = "ENSEMBL",toType = "ENTREZID",OrgDb = org.Hs.eg.db)
gene_ids <- gene_df$ENTREZID
head(gene_df)

#GO Enrichment Analysis
go_result <- enrichGO(gene=gene_ids,OrgDb=org.Hs.eg.db,keyType="ENTREZID",
                      ont="BP",pAdjustMethod="BH",qvalueCutoff=0.05)
head(go_result)

#Plot GO result
dotplot(go_result, showCategory=10)

#KEGG Pathway Analysis
kegg_result <- enrichKEGG(gene=gene_ids,organism="hsa",pvalueCutoff=0.05)
head(kegg_result)

#KEGG visualization
dotplot(kegg_result, showCategory=10)

#Survival Analysis
#Pick a gene
head(rownames(counts))
rownames(counts) <- gsub("\\..*","", rownames(counts))
tp53_id <- bitr("TP53",
                fromType = "SYMBOL",
                toType = "ENSEMBL",
                OrgDb = org.Hs.eg.db)$ENSEMBL
expr <- counts[tp53_id, ]
head(rownames(counts))

#Get survival data
clinical <- colData(data)
surv_data <- data.frame(expr=expr,time=clinical$days_to_death,status=clinical$vital_status)

#Clean data
surv_data <- na.omit(surv_data)

#Convert status
surv_data$status <- ifelse(surv_data$status == "Dead", 1, 0)

#Split groups
surv_data$group <- ifelse(surv_data$expr > median(surv_data$expr), "High","Low")

#Kaplan-Meier model
#install.packages(c("survival","survminer"))
library(survival)
library(survminer)
fit <- survfit(Surv(time, status)~group,data = surv_data)

#Plot survival curve
ggsurvplot(fit)