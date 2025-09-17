#!/bin/bash

# =============================================================================
# Cell Ranger Multi Processing Script for scRNA-seq Data
# =============================================================================
# Description: Processes FASTQ files using Cell Ranger Multi to generate 
#              feature-barcode matrices for single-cell RNA sequencing analysis
# Author: Hamid Khoshfekr Rudsari
# Email: khoshfekr1994@gmail.com
# Date: September 2025
# =============================================================================

# LSF Job Configuration
#BSUB -J cellranger_multi                    # Job name
#BSUB -n 28                                  # Number of CPU cores
#BSUB -M 350                                 # Maximum memory limit in GB
#BSUB -R "rusage[mem=16] span[hosts=1]"      # Memory per core (16GB) and keep on same host
#BSUB -W 200:00                              # Wall clock limit (200 hours for large datasets)
#BSUB -q long                                # Queue name (adjust for your HPC system)
#BSUB -o %J.cellranger.out                   # Standard output file
#BSUB -e %J.cellranger.err                   # Standard error file  
#BSUB -u your.email@institution.edu         # Email for notifications (CHANGE THIS)
#BSUB -N                                     # Send notification on job completion
#BSUB -B                                     # Send notification when job begins

# Set working directory (CHANGE THIS TO YOUR PROJECT PATH)
#BSUB -cwd /path/to/your/project/directory

# =============================================================================
# Configuration Variables (MODIFY THESE PATHS)
# =============================================================================

# Sample identifier (modify as needed)
SAMPLE_ID="S3"

# Configuration file path (ensure this file exists in your working directory)
CONFIG_FILE="multi_config.csv"

# =============================================================================
# Module Loading and Environment Setup
# =============================================================================

# Load Cell Ranger module (version may vary by system)
module load cellranger/9.0.0

# Print job information for debugging
echo "=========================================="
echo "Cell Ranger Multi Processing Job Started"
echo "=========================================="
echo "Job ID: $LSB_JOBID"
echo "Sample ID: $SAMPLE_ID" 
echo "Config File: $CONFIG_FILE"
echo "Working Directory: $(pwd)"
echo "Date: $(date)"
echo "Host: $(hostname)"
echo "=========================================="

# Verify configuration file exists
if [ ! -f "$CONFIG_FILE" ]; then
    echo "ERROR: Configuration file $CONFIG_FILE not found!"
    echo "Please ensure multi_config.csv exists in the working directory"
    exit 1
fi

# =============================================================================
# Cell Ranger Multi Execution
# =============================================================================

# Run Cell Ranger Multi
echo "Starting Cell Ranger Multi processing..."
echo "Command: cellranger multi --id $SAMPLE_ID --csv $CONFIG_FILE"

cellranger multi --id "$SAMPLE_ID" --csv "$CONFIG_FILE"

# Check if Cell Ranger completed successfully
if [ $? -eq 0 ]; then
    echo "=========================================="
    echo "Cell Ranger Multi completed successfully!"
    echo "Output directory: ${SAMPLE_ID}/"
    echo "Completion time: $(date)"
    echo "=========================================="
else
    echo "=========================================="
    echo "ERROR: Cell Ranger Multi failed!"
    echo "Check error logs for details"
    echo "Failure time: $(date)"
    echo "=========================================="
    exit 1
fi

# =============================================================================
# Post-processing Information
# =============================================================================

echo "Next steps:"
echo "1. Check the web summary report: ${SAMPLE_ID}/web_summary.html"
echo "2. Verify output files in: ${SAMPLE_ID}/sample_filtered_feature_bc_matrix/"
echo "3. Proceed with R/Seurat analysis using the generated matrix files"
echo "=========================================="
