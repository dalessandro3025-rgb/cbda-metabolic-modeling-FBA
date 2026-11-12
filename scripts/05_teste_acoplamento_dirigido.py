"""
Etapa 5: teste dirigido de acoplamento (hipótese do dreno de NADH).

Contexto
--------
A varredura das etapas 2-3 não encontrou acoplamento, mas cobriu apenas
deleções simples (exaustiva) e um subconjunto de duplas. O mecanismo usual
de acoplamento crescimento-produto é tornar a via do produto o ÚNICO dreno
disponível para um cofator em excesso (como ocorre com etanol e succinato).

A via do CBDA consome 4 NADH por molécula (HBD_c4/c6, TER_c4/c6), portanto
é candidata a dreno de NADH. A objeção levantada na etapa 3 foi que o
acoplamento por NADH exige anaerobiose, incompatível com a CBDAS (oxidase
O2-dependente). Mas há uma terceira via: manter O2 disponível para a CBDAS
e bloquear a reoxidação de NADH pelas rotas fermentativas E respiratórias.

Este script testa essa hipótese explicitamente, com knockouts combinados
que a varredura por pares não alcançaria.
"""

import sys
from itertools import combinations
from pathlib import Path

import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_model import build

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

BIO = "BIOMASS_SC5_notrace"
TARGET = "DM_cbda_c"
TOL = 1e-6

# Drenos alternativos de NADH em S. cerevisiae
GRUPOS = {
    "fermentacao_etanol": ["ADH1", "ADH2", "ADH3", "ADH4", "ADH5"],
    "glicerol": ["GPD1", "GPD2"],
    "nadh_desidrogenases_externas": ["NDE1", "NDE2"],
    "nadh_desidrogenase_interna": ["NDI1"],
    "aldeido_desidrogenase": ["ALD4", "ALD5", "ALD6"],
}


def ids_por_nome(modelo, nomes):
    achados = []
    for g in modelo.genes:
        if (g.name or "").upper() in {n.upper() for n in nomes}:
            achados.append(g.id)
    return achados


def avaliar(modelo):
    modelo.objective = BIO
    mu = modelo.slim_optimize()
    if mu is None or pd.isna(mu) or mu < TOL:
        return 0.0, 0.0, 0.0
    fva = flux_variability_analysis(modelo, reaction_list=[TARGET],
                                    fraction_of_optimum=0.999)
    return mu, float(fva.loc[TARGET, "minimum"]), float(fva.loc[TARGET, "maximum"])


def main():
    modelo = build()
    mu0, cmin0, cmax0 = avaliar(modelo)
    print(f"selvagem: mu={mu0:.4f} cbda_min={cmin0:.4f} cbda_max={cmax0:.4f}\n")

    mapa = {k: ids_por_nome(modelo, v) for k, v in GRUPOS.items()}
    for k, v in mapa.items():
        print(f"{k}: {len(v)} genes localizados")

    linhas = []
    chaves = list(mapa)
    # Todas as combinações de grupos funcionais (1 a 5 grupos simultâneos)
    for n in range(1, len(chaves) + 1):
        for combo in combinations(chaves, n):
            alvos = [g for k in combo for g in mapa[k]]
            if not alvos:
                continue
            with modelo as mm:
                for gid in alvos:
                    mm.genes.get_by_id(gid).knock_out()
                mu, cmin, cmax = avaliar(mm)
            linhas.append({
                "combinacao": " + ".join(combo),
                "n_genes": len(alvos),
                "mu": mu,
                "mu_rel": mu / mu0 if mu0 else 0.0,
                "cbda_min": cmin,
                "cbda_max": cmax,
                "acoplado": cmin > TOL,
            })

    df = pd.DataFrame(linhas).sort_values(["cbda_min", "mu"], ascending=False)
    df.to_csv(RES / "acoplamento_dirigido.csv", index=False)
    print("\n=== combinações de grupos funcionais ===")
    print(df.to_string(index=False))

    acop = df[df["acoplado"] & (df["mu_rel"] > 0.01)]
    print(f"\ncombinações com acoplamento e crescimento viável: {len(acop)}")

    # ------------------------------------------------------------------
    # Teste complementar: anaerobiose estrita (CBDAS sem O2) como controle
    # ------------------------------------------------------------------
    with modelo as mm:
        mm.reactions.EX_o2_e.lower_bound = 0
        mu, cmin, cmax = avaliar(mm)
    print(f"controle anaeróbio: mu={mu:.4f} cbda_min={cmin:.4f} cbda_max={cmax:.4f}")

    # ------------------------------------------------------------------
    # Limite teórico: o CBDA pode ser obrigatório sob ALGUMA restrição?
    # Força-se crescimento no máximo e minimiza-se CBDA sem nenhum knockout.
    # ------------------------------------------------------------------
    with modelo as mm:
        mm.reactions.get_by_id(BIO).lower_bound = 0.99 * mu0
        mm.objective = TARGET
        mm.objective_direction = "min"
        print(f"CBDA mínimo forçando mu >= 99% do maximo: {mm.slim_optimize():.6f}")


if __name__ == "__main__":
    main()
