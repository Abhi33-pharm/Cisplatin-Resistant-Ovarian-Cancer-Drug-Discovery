# Install the required packages

pak::pkg_install(c("dynamicTreeCut", "fastcluster", "flashClust"))


# Load the required library

library(tidyverse)
library(WGCNA)
library(data.table)
library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)
library(RColorBrewer)
library(pathview)
library(pheatmap)
library(EnhancedVolcano)
library(GEOquery)
library(AnnotationDbi)


# Enable multi-threading for WGCNA
enableWGCNAThreads()

# 2. Load and Format Expression & Trait Data

cat("Loading baseline expression and trait data...\n")

expr_data <- read.csv("07_baseline_expr_all.csv", 
                      row.names = 1, check.names = FALSE)

meta_data <- read.csv("08_baseline_meta_all.csv", 
                      row.names = 1, check.names = FALSE)

 

# In WGCNA, rows should be samples (cell lines) and columns MUST be genes
# Verify orientation:

if (ncol(expr_data) > nrow(expr_data)) {
  datExpr0 <- as.data.frame(expr_data) 
} else {
  datExpr0 <- as.data.frame(t(expr_data))
}

# Filter Invariant / Low-Variance Genes

N_GENES <- 5000
gene_vars <- apply(datExpr0, 2, var, na.rm = TRUE)
top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:min(N_GENES, length(gene_vars))]
datExpr <- datExpr0[, top_genes]

cat(sprintf("Retained top %d most variable genes across %d cell lines.\n", 
            ncol(datExpr), nrow(datExpr)))

# Check sample quality & missingness
gsg <- goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}

# Align metadata traits to datExpr rows

common_samples <- intersect(rownames(datExpr), rownames(meta_data))
datExpr <- datExpr[common_samples, ]
datTraits <- meta_data[common_samples, , drop = FALSE]


# Soft-Thresholding Power Selection
 
cat("Evaluating scale-free topology across soft-threshold powers...\n")

fig_dir=r"(C:\Users\naemi\Documents\Cis-platin\Cisplatin-Resistant-Ovarian-Cancer-Drug-Discovery\Results\Figures\plots)"

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

# Plot Scale Independence & Mean Connectivity
png(file.path(fig_dir, "06_soft_threshold_selection.png"), width = 2400, height = 1200, res = 300)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 2.5, 1))

# Scale-free topology fit index

r_sq <- -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2]
plot(sft$fitIndices[, 1], r_sq,
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit (signed R^2)",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], r_sq, labels = powers, col = "red", cex = 0.9)
abline(h = 0.85, col = "blue", lty = 2)

# Mean connectivity
plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)",
     ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, col = "darkgreen", cex = 0.9)

dev.off()



# Choosing power 
softPower <- sft$powerEstimate
if (is.na(softPower)) {
  softPower <- 12  
  cat(sprintf("Automatic power did not reach 0.85; defaulting to standard power = %d\n", softPower))
} else {
  cat(sprintf("Selected optimal soft-thresholding power: %d\n", softPower))
}


cor <- WGCNA::cor

# Run blockwiseModules
 
net <- blockwiseModules(
  datExpr,
  power = 12,                  
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,          
  deepSplit = 3,               
  mergeCutHeight = 0.10,      
  reassignThreshold = 0,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)
# Convert numeric labels to module colors
moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs

# Plot Gene Dendrogram with Module Colors
png(file.path(fig_dir, "07_gene_cluster_dendrogram.png"), width = 2600, height = 1400, res = 300)
plotDendroAndColors(
  net$dendrograms[[1]],
  moduleColors[net$blockGenes[[1]]],
  "Module Colors",
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  main = "Gene Clustering Dendrogram and Co-expression Modules"
)
dev.off()

 
# Module-Trait Correlation Analysis (GI50)
 
# calculate module eigengenes with color labels

MEs0 <- moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)

# Select numeric drug-response traits (GI50 and logGI50)

trait_cols <- intersect(c("GI50", "logGI50"), colnames(datTraits))
trait_subset <- as.data.frame(lapply(datTraits[, trait_cols, drop = FALSE], as.numeric))
rownames(trait_subset) <- rownames(datTraits)

# Calculate Pearson correlation and Student asymptotic p-values
moduleTraitCor <- cor(MEs, trait_subset, use = "pairwise.complete.obs")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# Format text for heatmap cells
textMatrix <- paste0(signif(moduleTraitCor, 2), "\n(p = ", signif(moduleTraitPvalue, 2), ")")
dim(textMatrix) <- dim(moduleTraitCor)

# Plot Module-Trait Relationships Heatmap
png(file.path(fig_dir, "08_module_trait_relationships.png"), width = 1800, height = 2400, res = 300)
par(mar = c(6, 9, 3, 3))
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(trait_subset),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.8,
  zlim = c(-1, 1),
  main = "Module-Trait Relationships (Cisplatin GI50)"
)
dev.off()


 
# Export Gene-Module Membership Table
 

out_table_dir=r"(C:\Users\naemi\Documents\Cis-platin\Cisplatin-Resistant-Ovarian-Cancer-Drug-Discovery\Results\Tables)"



gene_module_df <- data.frame(
  Gene = colnames(datExpr),
  Module = moduleColors,
  stringsAsFactors = FALSE
)

out_gene_module <- file.path(out_table_dir, "12_WGCNA_gene_module_assignments.csv")
write.csv(gene_module_df, out_gene_module, row.names = FALSE)

 


#  Direct Extraction of All Brown Module Genes
brown_genes <- colnames(datExpr)[moduleColors == "brown"]
cat(sprintf("Total genes in the brown module: %d\n", length(brown_genes)))


 
#  Calculating the Module Membership (MM / kME) and Gene Significance (GS)
 
# Module Membership: Correlation between individual gene expression and MEbrown eigengene
geneModuleMembership <- as.data.frame(cor(datExpr, MEs, use = "pairwise.complete.obs"))
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nrow(datExpr)))

# Gene Significance: Correlation between individual gene expression and logGI50 trait
geneTraitSignificance <- as.data.frame(cor(datExpr, trait_subset$logGI50, use = "pairwise.complete.obs"))
GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nrow(datExpr)))

names(geneTraitSignificance) <- "GS.logGI50"
names(GSPvalue) <- "p.GS.logGI50"

 
# Creating a Detailed Hub Gene Table for the Brown Module
 
brown_summary <- data.frame(
  Gene = brown_genes,
  Module = "brown",
  kME = geneModuleMembership[brown_genes, "MEbrown"],
  p_kME = MMPvalue[brown_genes, "MEbrown"],
  GS_logGI50 = geneTraitSignificance[brown_genes, "GS.logGI50"],
  p_GS = GSPvalue[brown_genes, "p.GS.logGI50"],
  stringsAsFactors = FALSE
)

# Rank genes by highest intra-modular connectivity (hub gene potential)
brown_summary <- brown_summary[order(-abs(brown_summary$kME)), ]

 
#  Save Extracted Gene Lists to CSV
 
# Full table with module connectivity metrics
out_brown_detailed <- file.path(out_table_dir, "13_brown_module_genes_detailed.csv")
write.csv(brown_summary, out_brown_detailed, row.names = FALSE)

# Simple list of gene symbols for functional enrichment (clusterProfiler / Reactome)
out_brown_list <- file.path(out_table_dir, "13_brown_module_gene_list.csv")
writeLines(brown_genes, out_brown_list)















