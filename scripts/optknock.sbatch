#!/bin/bash
#SBATCH --job-name=optknock_cbda
#SBATCH --partition=max50
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=48:00:00
#SBATCH --output=logs/optknock_%j.out
#SBATCH --error=logs/optknock_%j.err

# OptKnock para a via do CBDA em iMM904.
# Requer o ambiente descrito em docs/ambiente_cameo.md.

set -euo pipefail
mkdir -p logs results

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cameo_optknock

echo "host: $(hostname)"
echo "inicio: $(date -Is)"

python -u scripts/06_optknock.py

echo "fim: $(date -Is)"
