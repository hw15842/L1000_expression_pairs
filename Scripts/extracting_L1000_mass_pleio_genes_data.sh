#!/bin/bash

#SBATCH --job-name=l1000
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=20
#SBATCH --time=12:00:00
#SBATCH --mem=10000M
#SBATCH --partition=cpu


echo "Running on ${HOSTNAME}"
module add languages/r/4.0.3 
Rscript extracting_L1000_mass_pleio_genes_data.R /mnt/storage/scratch/hw15842/repo/L1000_expression_pairs/Data/