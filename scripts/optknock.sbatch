#!/bin/bash
#SBATCH --job-name=optknock_cbda
#SBATCH --partition=max50
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G
#SBATCH --output=logs/optknock_%j.out
#SBATCH --error=logs/optknock_%j.err

# Sem --time: herda o tempo padrao da particao (~3-4 dias na davinci).
# Dimensionamento medido: pico ~450 MB, single-thread (GLPK nao paraleliza).
# O custo e TEMPO, nao recurso. Nao pedir GPU (evita QOSMaxGRESPerUser).

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD}"
mkdir -p logs results

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate cameo_optknock

echo "host:   $(hostname)"
echo "job:    ${SLURM_JOB_ID:-local}"
echo "dir:    $PWD"
echo "python: $(python --version 2>&1)"
echo "inicio: $(date -Is)"
echo

python -u scripts/06_optknock.py 3

echo
echo "fim: $(date -Is)"
