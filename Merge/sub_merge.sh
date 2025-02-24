#!/bin/bash
#SBATCH --time=0-10:00:00
#SBATCH --account=def-salehlab-ab
#SBATCH --nodes=1
#SBATCH --output=Log_Output.out
#SBATCH --mem=600G
#SBATCH --cpus-per-task=30
#SBATCH --mail-user=fenosoa.randrianjatovo@aims.ac.rw
#SBATCH --mail-type=ALL

echo "Job started ..., Let's do it!!!"

module load StdEnv/2023
module spider arrow/18.1.0
module load arrow/18.1.0
module spider gdal
module spider proj
module spider geos
module spider udunits
module load  gcc/12.3  r/4.4.0


Rscript data_merge_analysis.R

echo "Job finished"
