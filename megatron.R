# =================================================================
# Environment Setup - D Drive because C is full of memes
# =================================================================
Sys.setenv(TMP = "D:/rtemp") # I dont have money to buy an ssd ngl
Sys.setenv(TEMP = "D:/rtemp")
Sys.setenv(TMPDIR = "D:/rtemp")
if (!dir.exists("D:/rtemp")) {
  dir.create("D:/rtemp") # Hope this works, had permission issues last week
  message("Created temp directory - hope this doesn't break later!")
} else {
  message("Temp directory already exists")
}

# For R packages 
if (!dir.exists("D:/R_libs")) {
  dir.create("D:/R_libs")
  message("Made new R library folder")
}
.libPaths("D:/R_libs") # Override default library location
message("\nCurrent library paths:", paste(.libPaths(), collapse = "\n"))

# =================================================================
# Package Installation - Following tutorial but confused about BiocManager
# =================================================================
if (!require("BiocManager", quietly = TRUE)) {
  message("\nInstalling BiocManager")
  install.packages("BiocManager") # Took 10 minutes last time
}

required_packages <- c("TCGAbiolinks", "DESeq2", "org.Hs.eg.db", "ggplot2", 
                       "pheatmap", "ggrepel", "survival", "survminer")

message("\nChecking installed packages...")
installed <- installed.packages()
missing <- required_packages[!required_packages %in% installed]
if(length(missing) > 0) {
  message("Missing packages: ", paste(missing, collapse = ", "))
  BiocManager::install(missing) 
  message("Installation attempt complete")
} else {
  message("All packages already installed - lfg")
}

# =================================================================
# Load Libraries 
# =================================================================
message("\nLoading libraries...")
library(TCGAbiolinks) # For TCGA data download
library(DESeq2)       # RNA-seq analysis
library(ggplot2)      # Plots
library(pheatmap)     # Heatmaps (v hot)
library(ggrepel)      # For labeling points
library(survival)     # Survival analysis
library(survminer)    # Survival plots
message("Libraries loaded")

# =================================================================
# Config
# =================================================================
cancer_type <- "TCGA-BRCA" # Chose BRCA because aunt had it
download_dir <- "D:/TCGA_data" 
output_dir <- "D:/TCGA_project/TCGA_results"


message("\nChecking output directory...")
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
  message("Created output directory at ", output_dir)
} else {
  message(output_dir, " already exists")
}

# =================================================================
# 1. Download & Prepare Expression Data 
# =================================================================
expression_data_path <- file.path(download_dir, paste0(cancer_type, "_expression.rds"))

if (!file.exists(expression_data_path)) {
  message("\nStarting data download... go get coffee!")
  query <- GDCquery(
    project = cancer_type,
    data.category = "Transcriptome Profiling",
    data.type = "Gene Expression Quantification", # Should be counts?
    workflow.type = "STAR - Counts" # Tutorial said use this
  )
  
  message("Downloading data... might take 1h+")
  GDCdownload(query, directory = download_dir) # 45 minutes later...
  
  message("Preparing data...")
  expression_data <- GDCprepare(query, directory = download_dir)
  
  message("Saving to RDS...")
  saveRDS(expression_data, expression_data_path)
  message("Data saved! File size:", file.size(expression_data_path)/1e6, "MB")
} else {
  message("\nLoading existing expression data from", expression_data_path)
  expression_data <- readRDS(expression_data_path)
  message("Data loaded! Dimensions:", dim(expression_data))
}

# Check what we actually got
message("\nSample metadata columns:")
print(colnames(colData(expression_data))[1:5]) # First 5 columns

message("\nFirst 3 sample IDs:")
print(head(rownames(colData(expression_data)), 3))

# =================================================================
# 2. Process Sample Types - Filtering is confusing
# =================================================================
message("\nProcessing sample types...")
sample_metadata <- as.data.frame(colData(expression_data))
count_matrix <- assay(expression_data, "unstranded")

message("Original sample types:")
print(table(sample_metadata$sample_type)) # See what we're working with

# Categorize sampl
sample_metadata <- sample_metadata %>%
  mutate(
    sample_type = case_when(
      grepl("Normal", sample_type) ~ "Normal",
      grepl("Tumor|Metastatic", sample_type) ~ "Tumor",
      TRUE ~ "Other" # For weird cases
    )
  ) %>% 
  filter(sample_type %in% c("Tumor", "Normal")) # Remove Others

message("\nFiltered sample counts:")
print(table(sample_metadata$sample_type))

# Match counts to filtered samples
count_matrix <- count_matrix[, rownames(sample_metadata)]
message("\nCount matrix dimensions after filtering:", dim(count_matrix))

# =================================================================
# 3. Differential Expression Analysis - DESeq2 
# =================================================================
message("\nCreating DESeq object...")
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ sample_type # Tumor vs Normal
)

message("Running DESeq - go take a nap...")
dds <- DESeq(dds) # 30+ minutes on my laptop

message("\nExtracting results...")
res <- results(dds, contrast = c("sample_type", "Tumor", "Normal"))

# Fix gene IDs 
message("\nOriginal gene IDs look like:", head(rownames(res)))
rownames(res) <- sub("\\..*", "", rownames(res))
message("Cleaned gene IDs:", head(rownames(res)))

# Map to gene symbols 
message("\nMapping ENSEMBL to symbols...")
res$symbol <- mapIds(org.Hs.eg.db,
                     keys = rownames(res),
                     column = "SYMBOL",
                     keytype = "ENSEMBL",
                     multiVals = "first")

message("Checking TP53 (should be present):")
tp53_row <- res[rownames(res) == "ENSG00000141510",]
print(paste("TP53 symbol:", tp53_row$symbol, 
            "log2FC:", round(tp53_row$log2FoldChange, 2)))

message("\nFirst 5 DEGs:")
print(head(res[order(res$padj),], 5))

# Save results - Excel can't handle big files
write.csv(as.data.frame(res), file.path(output_dir, "DEG_results.csv"))
message("\nSaved DEG results to CSV")

# =================================================================
# 4. Visualization - its shit but ill work on it 
# =================================================================
# Volcano plot 
message("\nPreparing volcano plot...")
res_clean <- res[!is.na(res$log2FoldChange) & !is.na(res$pvalue),]
sig_res <- res_clean[res_clean$padj < 0.05 & abs(res_clean$log2FoldChange) > 1,]

message("Significant genes:", nrow(sig_res))
if(nrow(sig_res) == 0) stop("No significant genes")

top_genes <- head(sig_res[order(sig_res$padj),], 20)
message("Top gene:", top_genes$symbol[1])

volcano_plot <- ggplot(sig_res, aes(x = log2FoldChange, y = -log10(pvalue))) +
  geom_point(aes(color = "Significant"), alpha = 0.7) +
  geom_text_repel(data = top_genes, aes(label = symbol), 
                  max.overlaps = 20, size = 3) + # Smaller text
  scale_color_manual(values = "red") +
  labs(title = "Tumor vs Normal DEGs",
       subtitle = paste(nrow(sig_res), "significant genes"),
       x = "log2 Fold Change",
       y = "-log10 p-value") +
  theme_minimal()

ggsave(file.path(output_dir, "volcano_plot.png"), volcano_plot, 
       width = 10, height = 7, dpi = 300)

# Heatmap - gene names never fit
message("\nCreating heatmap...")
top_genes <- rownames(res)[head(order(res$padj), 50)]
normalized_counts <- counts(dds, normalized = TRUE)

# Fix row names for genes
rownames(normalized_counts) <- sub("\\..*", "", rownames(normalized_counts))
heatmap_data <- normalized_counts[top_genes,]
rownames(heatmap_data) <- res[top_genes, "symbol"] # Use symbols

annotation_df <- data.frame(SampleType = sample_metadata$sample_type)
rownames(annotation_df) <- colnames(heatmap_data)

message("Heatmap dimensions:", dim(heatmap_data))
if(nrow(heatmap_data) == 0) stop("No genes for heatmap")

pheatmap(log2(heatmap_data + 1),
         annotation_col = annotation_df,
         show_rownames = TRUE, # Risk overlapping text
         fontsize_row = 7,     # Tiny but readable?
         main = "Top 50 DEGs Heatmap",
         filename = file.path(output_dir, "heatmap.png"))

message("Saved heatmap - hope labels are visible!")

# PCA plot 
message("\nRunning PCA...")
vsd <- vst(dds, blind = FALSE)
pca_data <- plotPCA(vsd, intgroup = "sample_type", returnData = TRUE)

message("PCA variance explained:")
print(attr(pca_data, "percentVar"))

pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = sample_type)) +
  geom_point(size = 3) +
  stat_ellipse(level = 0.95) + 
  labs(title = "PCA Plot",
       x = paste0("PC1 (", round(100*attr(pca_data, "percentVar")[1],1), "%)"),
       y = paste0("PC2 (", round(100*attr(pca_data, "percentVar")[2],1), "%)")) +
  theme_minimal()

ggsave(file.path(output_dir, "pca_plot.png"), pca_plot, 
       width = 9, height = 6)

# =================================================================
# 5. Survival Analysis - Not sure if doing this right
# =================================================================
message("\nAttempting survival analysis...")
clinical_data <- GDCquery_clinic(project = cancer_type, type = "clinical")

# Merge with expression data
common_patients <- intersect(substr(rownames(sample_metadata), 1, 12), 
                             clinical_data$bcr_patient_barcode)
clinical_filtered <- clinical_data[clinical_data$bcr_patient_barcode %in% common_patients,]

# Get TP53 expression
tp53_counts <- log2(normalized_counts["ENSG00000141510",] + 1)
tp53_df <- data.frame(
  patient_barcode = substr(names(tp53_counts), 1, 12),
  tp53_expression = tp53_counts
)

survival_df <- merge(clinical_filtered, tp53_df, by = "bcr_patient_barcode")

message("\nPatients in survival analysis:", nrow(survival_df))

# Kaplan-Meier plot 
surv_obj <- Surv(time = survival_df$days_to_last_follow_up,
                 event = survival_df$vital_status == "Dead")

surv_fit <- survfit(surv_obj ~ (tp53_expression > median(tp53_expression)),
                    data = survival_df)

surv_plot <- ggsurvplot(surv_fit, 
                        data = survival_df,
                        pval = TRUE,
                        title = "TP53 Expression Survival Analysis",
                        xlab = "Days",
                        risk.table = TRUE)

ggsave(file.path(output_dir, "survival_plot.png"), surv_plot$plot,
       width = 10, height = 8)

# =================================================================
# 6. ML Prep - No idea what algorithm to use yet
# =================================================================
message("\nPreparing ML data...")
sig_genes <- rownames(res)[res$padj < 0.05 & abs(res$log2FoldChange) > 1]
top_100 <- head(sig_genes[order(res[sig_genes, "padj"])], 100)

ml_data <- data.frame(t(normalized_counts[top_100, ]))
ml_data$sample_type <- sample_metadata$sample_type

message("ML data dimensions:", dim(ml_data))
message("Class balance:")
print(table(ml_data$sample_type))

write.csv(ml_data, file.path(output_dir, "ml_data.csv"), row.names = FALSE)
message("Saved ML data")

# =================================================================
# 7. Session Info - For reproducibility
# =================================================================
sink(file.path(output_dir, "session_info.txt"))
message("\nSession info for reproducibility:")
print(sessionInfo())
sink()
message("Analysis complete! Time:", Sys.time())

getwd()
setwd("D:/BIN_TCGA_BCA")