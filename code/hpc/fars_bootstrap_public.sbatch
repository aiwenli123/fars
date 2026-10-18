#!/bin/bash
#SBATCH --job-name=fars-bootstrap
#SBATCH --account=YOUR_PROJECT_ACCOUNT
#SBATCH --partition=genoa-std-mem
#SBATCH --array=1-100
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=12:00:00
#SBATCH --output=logs/fars_%A_%a.out
#SBATCH --error=logs/fars_%A_%a.err

set -e

# Load Wharton software and R.
source /vast/projects/wharton/software/bin/wharton_parcc.sh
module load rstudio

# Run the simulation from the project directory.
cd /path/to/fars/project
Rscript --vanilla 09a_simulation_full.R