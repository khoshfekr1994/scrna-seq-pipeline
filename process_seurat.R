# =============================================================================
# Single-cell RNA-seq Data Processing with Seurat
# =============================================================================
# Description: Complete pipeline for processing Cell Ranger output through
#              quality control, normalization, dimensionality reduction, 
#              and clustering analysis
# Author: Hamid Khoshfekr Rudsari
# Email: khoshfekr1994@gmail.com
# Date: September 2025
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
})

# =============================================================================
# Configuration Parameters
# =============================================================================

# File paths (MODIFY THESE PATHS)
input_dir <- "path/to/cellranger/output/sample_filtered_feature_bc_matrix/"
output_dir <- "path/to/output/directory/"

# Sample information
sample_name <- "S1"  # Change to your sample name
project_name <- "scRNA_project"  # Change to your project name

# Quality control thresholds
min_cells_per_gene <- 3      # Include genes detected in at least this many cells
min_features_per_cell <- 200 # Include cells with at least this many genes
max_features_per_cell <- 2500 # Maximum genes per cell (filter potential doublets)
max_mt_percent <- 20         # Maximum mitochondrial gene percentage

# Analysis parameters
pca_dims <- 1:30            # PCA dimensions to use for downstream analysis
clustering_resolution <- 0.5 # Clustering resolution (higher = more clusters)

# =============================================================================
# Helper Functions
# =============================================================================

# Function to create output directory if it doesn't exist
create_output_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE)
    cat("Created output directory:", path, "\n")
  }
}

# Function to save Seurat object with informative message
save_seurat_object <- function(seurat_obj, filename, description) {
  filepath <- file.path(output_dir, filename)
  saveRDS(seurat_obj, filepath)
  cat("Saved", description, "to:", filepath, "\n")
  cat("Object contains", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n\n")
}

# =============================================================================
# Data Loading and Initial Object Creation
# =============================================================================

cat("=================================================\n")
cat("Starting scRNA-seq data processing pipeline\n")
cat("=================================================\n")
cat("Sample:", sample_name, "\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# Create output directory
create_output_dir(output_dir)

# Load 10X data
cat("Loading 10X data...\n")
if (!dir.exists(input_dir)) {
  stop("ERROR: Input directory does not exist: ", input_dir)
}

data <- Read10X(data.dir = input_dir)
cat("Loaded expression matrix with", nrow(data), "genes and", ncol(data), "barcodes\n\n")

# Create Seurat object
cat("Creating initial Seurat object...\n")
seurat_obj <- CreateSeuratObject(
  counts = data, 
  project = project_name,
  min.cells = min_cells_per_gene,
  min.features = min_features_per_cell
)

cat("Initial object created with", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes\n\n")

# =============================================================================
# Quality Control Metrics
# =============================================================================

cat("Calculating quality control metrics...\n")

# Calculate mitochondrial gene percentage
# Note: Use "^MT-" for human, "^mt-" for mouse
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

# Calculate ribosomal gene percentage (optional)
seurat_obj[["percent.rb"]] <- PercentageFeatureSet(seurat_obj, pattern = "^Rp[sl]")

# Add number of genes and UMIs to metadata  
seurat_obj[["nCount_RNA"]] <- seurat_obj$nCount_RNA
seurat_obj[["nFeature_RNA"]] <- seurat_obj$nFeature_RNA

# Print QC summary statistics
cat("Quality Control Summary:\n")
cat("Mean genes per cell:", round(mean(seurat_obj$nFeature_RNA), 2), "\n")
cat("Mean UMIs per cell:", round(mean(seurat_obj$nCount_RNA), 2), "\n") 
cat("Mean mitochondrial %:", round(mean(seurat_obj$percent.mt), 2), "\n\n")

# =============================================================================
# Cell Filtering
# =============================================================================

cat("Applying quality control filters...\n")
cat("Filters applied:\n")
cat("- Minimum features per cell:", min_features_per_cell, "\n")
cat("- Maximum features per cell:", max_features_per_cell, "\n") 
cat("- Maximum mitochondrial %:", max_mt_percent, "\n\n")

# Apply filters
seurat_obj <- subset(
  seurat_obj, 
  subset = nFeature_RNA > min_features_per_cell & 
           nFeature_RNA < max_features_per_cell & 
           percent.mt < max_mt_percent
)

cat("After filtering:", ncol(seurat_obj), "cells and", nrow(seurat_obj), "genes retained\n\n")

# Save filtered object
save_seurat_object(seurat_obj, paste0(sample_name, "_filtered.rds"), "filtered Seurat object")

# =============================================================================
# Normalization and Scaling  
# =============================================================================

cat("Performing normalization with SCTransform...\n")

# SCTransform normalization (replaces NormalizeData, FindVariableFeatures, and ScaleData)
seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)

cat("SCTransform normalization completed\n\n")

# =============================================================================
# Dimensionality Reduction
# =============================================================================

cat("Performing dimensionality reduction...\n")

# Principal Component Analysis
cat("Running PCA...\n")
seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

# UMAP for visualization
cat("Running UMAP...\n")
seurat_obj <- RunUMAP(seurat_obj, dims = pca_dims, verbose = FALSE)

# Find neighbors for clustering
cat("Finding neighbors...\n")
seurat_obj <- FindNeighbors(seurat_obj, dims = pca_dims, verbose = FALSE)

cat("Dimensionality reduction completed\n\n")

# Save processed object (ideal checkpoint before clustering)
save_seurat_object(seurat_obj, paste0(sample_name, "_processed_ready_for_clustering.rds"), 
                   "processed Seurat object (ready for clustering)")

# =============================================================================
# Clustering
# =============================================================================

cat("Performing clustering...\n")
cat("Clustering resolution:", clustering_resolution, "\n")

# Find clusters
seurat_obj <- FindClusters(seurat_obj, resolution = clustering_resolution, verbose = FALSE)

# Get number of clusters
n_clusters <- length(unique(Idents(seurat_obj)))
cat("Found", n_clusters, "clusters\n\n")

# Save final clustered object
save_seurat_object(seurat_obj, paste0(sample_name, "_final_clustered.rds"), 
                   "final clustered Seurat object")

# =============================================================================
# Generate Visualizations
# =============================================================================

cat("Generating visualizations...\n")

# UMAP plot with clusters
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE, pt.size = 0.5) +
  ggtitle(paste("UMAP Clustering -", sample_name)) +
  theme(plot.title = element_text(hjust = 0.5))

# Feature plots for QC metrics
p2 <- FeaturePlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                  ncol = 3, pt.size = 0.3)

# Save plots
ggsave(file.path(output_dir, paste0(sample_name, "_umap_clusters.png")), 
       plot = p1, width = 8, height = 6, dpi = 300)
ggsave(file.path(output_dir, paste0(sample_name, "_qc_features.png")), 
       plot = p2, width = 15, height = 5, dpi = 300)

cat("Plots saved to output directory\n\n")

# =============================================================================
# Analysis Summary
# =============================================================================

cat("=================================================\n")
cat("Analysis completed successfully!\n")
cat("=================================================\n")
cat("Final Summary:\n")
cat("- Sample:", sample_name, "\n")
cat("- Total cells:", ncol(seurat_obj), "\n")
cat("- Total genes:", nrow(seurat_obj), "\n") 
cat("- Number of clusters:", n_clusters, "\n")
cat("- Output files saved to:", output_dir, "\n")
cat("\nGenerated files:\n")
cat("1.", paste0(sample_name, "_filtered.rds"), "- After QC filtering\n")
cat("2.", paste0(sample_name, "_processed_ready_for_clustering.rds"), "- Normalized and reduced\n")
cat("3.", paste0(sample_name, "_final_clustered.rds"), "- Final clustered object\n") 
cat("4.", paste0(sample_name, "_umap_clusters.png"), "- UMAP cluster visualization\n")
cat("5.", paste0(sample_name, "_qc_features.png"), "- QC metrics visualization\n")
cat("=================================================\n")

# Optional: Display plots if running interactively
if (interactive()) {
  print(p1)
  print(p2)
}
