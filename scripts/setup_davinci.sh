#!/bin/bash
# Cria o ambiente isolado para o OptKnock. Rodar UMA vez, da raiz do repo:
#     bash scripts/setup_davinci.sh
set -euo pipefail

ENV_NAME=cameo_optknock
source "$(conda info --base)/etc/profile.d/conda.sh"

if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    echo "ambiente ${ENV_NAME} ja existe, pulando a criacao"
else
    echo ">> criando ${ENV_NAME} (Python 3.11)"
    conda create -n "${ENV_NAME}" python=3.11 -y
fi

conda activate "${ENV_NAME}"

echo ">> instalando dependencias"
pip install --quiet "numpy<2" "cobra==0.26.3" cameo pip

echo ">> verificando"
python - <<'PY'
import sys
from pathlib import Path
sys.path.insert(0, str(Path("scripts").resolve()))
import cameo_compat
import cobra, cameo, numpy
from cameo.strain_design.deterministic.linear_programming import OptKnock
print(f"   numpy {numpy.__version__} | cobra {cobra.__version__} | cameo {cameo.__version__}")
print("   OptKnock importado com sucesso")
PY

echo
echo "pronto. Para submeter:  sbatch scripts/optknock.sbatch"
