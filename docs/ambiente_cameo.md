# Ambiente para OptKnock/OptGene (cameo)

O cameo (última versão 0.11.2) não é mais mantido e **não funciona** com o
ambiente principal do projeto (`fba_cbda`, COBRApy 0.31.1). As incompatibilidades
encontradas, em ordem:

| Erro | Causa |
|---|---|
| `ImportError: cannot import name 'trapz' from 'numpy'` | `numpy.trapz` removido no numpy 2.0 |
| `ImportError: cannot import name 'Mapping' from 'collections'` | `grako` (dependência de `gnomic`) é código Python 2; `collections.Mapping` removido no Python 3.10 |
| `ModuleNotFoundError: No module named 'pip'` / `pip has no attribute 'get_installed_distributions'` | API removida no pip 10 |
| `module 'sympy' has no attribute 'singleton'` | reorganização do sympy |
| `cannot import name 'find_gene_knockout_reactions'` | função removida no COBRApy ≥ 0.26 |
| `pandas.core.common has no attribute 'in_ipython_frontend'` | API movida no pandas |

Por isso o OptKnock roda em **ambiente separado**, sem tocar no `fba_cbda`.

## Criação do ambiente

```bash
# com conda (equivalente ao fba_cbda, mas isolado)
conda create -n cameo_optknock python=3.11 -y
conda activate cameo_optknock
pip install "numpy<2" "cobra==0.26.3" cameo pip
```

Os *shims* de compatibilidade estão em `scripts/cameo_compat.py` e precisam ser
importados **antes** de qualquer import do cameo — o `06_optknock.py` já faz isso.

## Execução

```bash
conda activate cameo_optknock
python scripts/06_optknock.py
```

## Custo computacional

O OptKnock é um problema de programação linear inteira mista (MILP) de dois
níveis. Com o solver GLPK (open source) e um espaço de busca de 141 reações
deletáveis, a execução não terminou em 30 minutos de CPU em máquina local.
O próprio cameo emite o aviso:

> `You are trying to run OptKnock with GLPK. This might not end well.`

Recomendações, em ordem de preferência:

1. **Rodar na davinci via Slurm** (ver `scripts/optknock.sbatch`), com parede de
   tempo generosa. É a opção mais prática.
2. **Usar um solver comercial** (Gurobi ou CPLEX, ambos com licença acadêmica
   gratuita). A diferença de desempenho em MILP é de ordens de magnitude.
   Depois de instalar: `modelo.solver = "gurobi"` antes de construir o OptKnock.
3. **Reduzir o espaço de busca**, restringindo `SUBSISTEMAS` em
   `scripts/06_optknock.py` ao metabolismo central estrito.
