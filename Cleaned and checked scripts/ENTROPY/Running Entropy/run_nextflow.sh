#!/bin/bash
#SBATCH --account=def-mcfarlas
#SBATCH --time=72:00:00  # Shorten it to 3 days
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --job-name=entropy_k2-9
#SBATCH --mail-user=lebla179@yorku.ca
#SBATCH --mail-type=ALL

# ------------------------------------------------------------
# Shell script for executing the main.nf script found in "running-entropy-mainnf.txt"
# Author: Jenna LeBlanc
# Created: August 2026
# Last Updated: March 10, 2026

# Note for users:
# Make sure main.nf and nextflow.config are created and in the same directory in order to run the pipeline
# ------------------------------------------------------------

# Load required modules
module load nextflow

# Run the pipeline
cd /home/lebla179/scratch/entropy_k_2-9/fixed_MPGL
nextflow run main.nf