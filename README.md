# scRNA-seq Processing Pipeline

A comprehensive pipeline for processing single-cell RNA sequencing data from FASTQ files to analysis-ready Seurat objects using Cell Ranger and R.

## Overview

This pipeline converts raw FASTQ files from single-cell RNA sequencing experiments into processed Seurat objects stored as RDS files. The workflow includes:

1. **Cell Ranger Multi**: Raw FASTQ processing and feature-barcode matrix generation
2. **R Processing**: Quality control, normalization, dimensionality reduction, and clustering using Seurat

## Prerequisites

### HPC Environment
- LSF batch system (for job submission)
- Cell Ranger 9.0.0 or compatible version
- Sufficient computational resources (recommended: 28+ cores, 350GB+ memory)

### R Dependencies
```r
# Install required packages
install.packages(c("Seurat", "dplyr", "ggplot2"))
```

### Required Files
- Raw FASTQ files from your scRNA-seq experiment
- `multi_config.csv` configuration file for Cell Ranger Multi
- Reference genome (human/mouse) compatible with Cell Ranger

## File Structure

```
project_directory/
├── cellranger_multi.sh          # Cell Ranger batch script
├── process_seurat.R             # R processing script  
├── multi_config.csv             # Cell Ranger configuration
├── fastq_files/                 # Raw FASTQ files
└── output/                      # Generated outputs
    ├── S3/                      # Cell Ranger output
    │   └── sample_filtered_feature_bc_matrix/
    ├── S1_filtered.rds          # QC filtered Seurat object
    ├── S1_processed_ready_for_clustering.rds  # Normalized object
    └── S1_final_clustered.rds   # Final clustered object
```

## Usage

### Step 1: Cell Ranger Multi Processing

1. **Prepare your configuration file** (`multi_config.csv`):
   ```csv
   [gene-expression]
   reference,/path/to/reference/genome
   fastqs,/path/to/fastq/files
   sample,sample_name
   ```

2. **Submit the Cell Ranger job**:
   ```bash
   bsub < cellranger_multi.sh
   ```

### Step 2: R Processing and Analysis

After Cell Ranger completes successfully, run the R processing script:

```bash
Rscript process_seurat.R
```

## Scripts Description

### `cellranger_multi.sh`
HPC batch script that:
- Allocates computational resources (28 cores, 350GB memory)
- Sets up appropriate queue and time limits
- Loads Cell Ranger module
- Runs `cellranger multi` with your configuration

**Key Parameters:**
- **Cores**: 28 (adjust based on your system)
- **Memory**: 350GB total, 16GB per core
- **Time Limit**: 200 hours (for large datasets)
- **Queue**: `long` (modify according to your HPC setup)

### `process_seurat.R`
R script that processes Cell Ranger output through multiple stages:

#### Quality Control & Filtering
- Loads 10X data using `Read10X()`
- Creates initial Seurat object with basic filtering
- Calculates mitochondrial gene percentages
- Filters cells based on:
  - Feature count: 200-2500 genes per cell
  - Mitochondrial content: <20%

#### Normalization & Dimensionality Reduction
- **SCTransform**: Modern normalization method
- **PCA**: Principal component analysis
- **UMAP**: Uniform manifold approximation and projection
- **Clustering**: Graph-based clustering with resolution 0.5

#### Checkpoints
The script creates three RDS files for different analysis stages:
1. `*_filtered.rds`: After QC filtering
2. `*_processed_ready_for_clustering.rds`: After normalization and dimension reduction
3. `*_final_clustered.rds`: Final clustered object

## Output Files

### Cell Ranger Output
- **Feature-barcode matrix**: Gene expression counts per cell
- **Metrics**: Summary statistics and QC metrics
- **Web summary**: HTML report with analysis overview

### Seurat RDS Files
- **Filtered object**: Quality-controlled cells and genes
- **Processed object**: Normalized data ready for clustering analysis
- **Final object**: Complete analysis with cluster assignments

### Visualizations
- **UMAP plot**: Cell clusters in reduced dimensional space
- **Feature plots**: QC metrics (gene count, UMI count, mitochondrial percentage)

## Customization

### Cell Ranger Parameters
Modify the batch script parameters based on your:
- **Dataset size**: Adjust memory and core allocation
- **HPC system**: Change queue names, time limits
- **File paths**: Update working directories and output paths

### R Analysis Parameters
Key parameters you might want to adjust:

```r
# Cell filtering thresholds
min.cells = 3           # Minimum cells per gene
min.features = 200      # Minimum genes per cell
max.features = 2500     # Maximum genes per cell (filter doublets)
percent.mt = 20         # Maximum mitochondrial percentage

# Analysis parameters  
dims = 1:30            # Number of PCA dimensions
resolution = 0.5       # Clustering resolution
```

## Troubleshooting

### Common Issues

1. **Memory errors**: Increase memory allocation in batch script
2. **Time limit exceeded**: Extend wall clock time for large datasets
3. **Missing reference**: Ensure Cell Ranger reference genome is properly installed
4. **R package errors**: Verify all required packages are installed

### Resource Estimation

| Dataset Size | Estimated Time | Memory Required |
|--------------|----------------|-----------------|
| ~5K cells    | 2-4 hours      | 150GB          |
| ~10K cells   | 4-8 hours      | 250GB          |
| ~20K+ cells  | 8+ hours       | 350GB+         |


## Author

Hamid Khoshfekr Rudsari, PhD
contact: 
khoshfekr1994@gmail.com
hkhoshfekr@mdanderson.org

## License

This project is licensed under the MIT License - see the LICENSE file for details.
