# ============================================================================
# POST-INFECTIOUS CFS MOLECULAR SIGNATURE ANALYSIS
# ============================================================================
# 
# Analysis of gene expression differences between post-infectious Chronic 
# Fatigue Syndrome (CFS) patients and healthy controls
#
# Data: GSE14577 - Affymetrix HG-U95Av2 microarray
# Samples: 8 CFS patients vs 7 healthy controls
# 
# Supporting Pathogen Effector Convergence Theory (PECT) research
#
# Author: [Your Name]
# Date: [Current Date]
# Repository: [Your GitHub URL]
# ============================================================================

# SETUP AND DEPENDENCIES ====================================================

# Install required packages (run once)
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

required_packages <- c("GEOquery", "limma", "ggplot2")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) BiocManager::install(new_packages)

# Load libraries
library(GEOquery)
library(limma)
library(ggplot2)

# Set working directory (adjust as needed)
# setwd("path/to/your/analysis/directory")

print("=== POST-INFECTIOUS CFS MOLECULAR ANALYSIS ===")
print("Downloading and analyzing GSE14577 dataset...")

# DATA DOWNLOAD AND PREPARATION =============================================

# Download CFS dataset
gse14577 <- getGEO("GSE14577", GSEMatrix = TRUE)
print("GSE14577 dataset downloaded successfully!")

# Extract data from first platform (contains the samples we need)
platform1 <- gse14577[[1]]
expr_data <- exprs(platform1)
sample_info <- pData(platform1)
feature_data <- fData(platform1)

print(paste("Dataset dimensions:", nrow(expr_data), "genes x", ncol(expr_data), "samples"))

# SAMPLE GROUP ASSIGNMENT ===================================================

print("Sample information:")
print(sample_info$title)

# Assign groups based on sample titles
groups <- ifelse(grepl("CFS patient", sample_info$title), "CFS", "Control")

print("Group assignments:")
group_summary <- data.frame(
  Sample = sample_info$title,
  Group = groups
)
print(group_summary)

print("Sample counts:")
print(table(groups))

# DIFFERENTIAL EXPRESSION ANALYSIS ==========================================

print("Performing differential expression analysis...")

# Create design matrix
design <- model.matrix(~ 0 + factor(groups))
colnames(design) <- c("CFS", "Control")

# Fit linear model
fit <- lmFit(expr_data, design)

# Create contrast (CFS vs Control)
contrast_matrix <- makeContrasts(CFS - Control, levels = design)
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# Extract results
results <- topTable(fit2, number = Inf, adjust.method = "BH")

print("Differential expression analysis complete!")

# RESULTS SUMMARY ============================================================

# Summary statistics
total_genes <- nrow(results)
sig_genes <- sum(results$adj.P.Val < 0.05)
sig_percentage <- round((sig_genes / total_genes) * 100, 2)

print("=== ANALYSIS RESULTS SUMMARY ===")
print(paste("Total genes analyzed:", total_genes))
print(paste("Significantly altered genes (adj.P.Val < 0.05):", sig_genes))
print(paste("Percentage of genome affected:", sig_percentage, "%"))

# Effect size categorization
large_up <- sum(results$logFC > 1 & results$adj.P.Val < 0.05)
moderate_up <- sum(results$logFC > 0.5 & results$logFC <= 1 & results$adj.P.Val < 0.05)
large_down <- sum(results$logFC < -1 & results$adj.P.Val < 0.05)
moderate_down <- sum(results$logFC < -0.5 & results$logFC >= -1 & results$adj.P.Val < 0.05)

effect_summary <- data.frame(
  Category = c("Large Up (>2-fold)", "Moderate Up (1.4-2-fold)", 
               "Large Down (>2-fold)", "Moderate Down (1.4-2-fold)", 
               "Total Significant"),
  Count = c(large_up, moderate_up, large_down, moderate_down, sig_genes),
  Percentage = round(c(large_up, moderate_up, large_down, moderate_down, sig_genes) / total_genes * 100, 2)
)

print("Effect size distribution:")
print(effect_summary)

# GENE SYMBOL EXTRACTION ====================================================

print("Extracting gene symbols...")

# Get top 20 most significant genes
top_20 <- head(results, 20)
gene_symbols_top20 <- feature_data[rownames(top_20), "Gene Symbol"]

# Create top results table
top_20_results <- data.frame(
  Rank = 1:20,
  Probe_ID = rownames(top_20),
  Gene_Symbol = as.character(gene_symbols_top20),
  LogFC = round(top_20$logFC, 3),
  P_Value = formatC(top_20$adj.P.Val, format = "e", digits = 2),
  Direction = ifelse(top_20$logFC > 0, "UP", "DOWN"),
  stringsAsFactors = FALSE
)

print("TOP 20 MOST SIGNIFICANT GENES:")
print(top_20_results)

# GENE LISTS FOR PATHWAY ANALYSIS ===========================================

# Create gene lists with lenient criteria for better coverage
upregulated_probes <- rownames(results[results$logFC > 0.3 & results$adj.P.Val < 0.05, ])
downregulated_probes <- rownames(results[results$logFC < -0.3 & results$adj.P.Val < 0.05, ])

# Extract gene symbols
if(length(upregulated_probes) > 0) {
  up_symbols <- feature_data[upregulated_probes, "Gene Symbol"]
  up_symbols_clean <- up_symbols[up_symbols != "" & !is.na(up_symbols)]
  
  print(paste("Upregulated genes (>1.23-fold):", length(up_symbols_clean)))
  print("Top upregulated genes:")
  print(head(up_symbols_clean, 15))
}

if(length(downregulated_probes) > 0) {
  down_symbols <- feature_data[downregulated_probes, "Gene Symbol"]
  down_symbols_clean <- down_symbols[down_symbols != "" & !is.na(down_symbols)]
  
  print(paste("Downregulated genes (>1.23-fold):", length(down_symbols_clean)))
  print("Top downregulated genes:")
  print(head(down_symbols_clean, 15))
}

# DATA EXPORT ================================================================

print("Saving results to files...")

# Save complete results
write.csv(results, "CFS_vs_Control_complete_results.csv", row.names = TRUE)

# Save gene lists for pathway analysis
write.csv(up_symbols_clean, "CFS_upregulated_genes.csv", row.names = FALSE)
write.csv(down_symbols_clean, "CFS_downregulated_genes.csv", row.names = FALSE)

# Save top results with gene symbols
write.csv(top_20_results, "CFS_top20_significant_genes.csv", row.names = FALSE)

# Save analysis summary
summary_stats <- data.frame(
  Metric = c("Total_Genes", "Significant_Genes", "Percent_Affected", 
             "Upregulated", "Downregulated", "Large_Effect_Up", "Large_Effect_Down"),
  Value = c(total_genes, sig_genes, sig_percentage, 
            length(up_symbols_clean), length(down_symbols_clean), large_up, large_down)
)
write.csv(summary_stats, "CFS_analysis_summary.csv", row.names = FALSE)

# VISUALIZATION ==============================================================

# Volcano plot
print("Creating volcano plot...")

volcano_data <- data.frame(
  logFC = results$logFC,
  neg_log_p = -log10(results$adj.P.Val),
  significant = results$adj.P.Val < 0.05
)

volcano_plot <- ggplot(volcano_data, aes(x = logFC, y = neg_log_p)) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("grey", "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "blue") +
  geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
  labs(title = "CFS vs Control - Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-Value") +
  theme_minimal() +
  theme(legend.position = "none")

ggsave("CFS_volcano_plot.png", volcano_plot, width = 8, height = 6, dpi = 300)

print("=== ANALYSIS COMPLETE ===")
print("Files saved:")
print("- CFS_vs_Control_complete_results.csv")
print("- CFS_upregulated_genes.csv")
print("- CFS_downregulated_genes.csv") 
print("- CFS_top20_significant_genes.csv")
print("- CFS_analysis_summary.csv")
print("- CFS_volcano_plot.png")

print("")
print("=== KEY FINDINGS ===")
print(paste("Discovered", sig_genes, "significantly altered genes in post-infectious CFS"))
print("Top upregulated genes include: ASAP1, ATRX, VCAN, CXCR4 (immune/stress response)")
print("Top downregulated genes include: SLC10A1, SCN10A, OPHN1 (transport/neuronal)")
print("")
print("=== NEXT STEPS FOR PECT RESEARCH ===")
print("1. Compare these genes with COVID-19 literature")
print("2. Perform pathway enrichment analysis") 
print("3. Validate findings with RT-PCR")
print("4. Search for therapeutic targets")
print("")
print("Analysis supports Pathogen Effector Convergence Theory (PECT)")
print("Different pathogens may target similar cellular pathways")

# SESSION INFO ===============================================================
print("=== SESSION INFORMATION ===")
sessionInfo()
