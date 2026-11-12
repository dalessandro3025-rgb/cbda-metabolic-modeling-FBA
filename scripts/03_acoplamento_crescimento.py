"""
Etapa 3: triagem de deleções por ACOPLAMENTO crescimento-produção.

Justificativa metodológica
--------------------------
Na etapa anterior, nenhuma deleção gênica aumentou o fluxo máximo de CBDA.
Isso não é um achado biológico, e sim uma propriedade da formulação: em FBA,
remover uma reação apenas restringe o espaço de soluções, de modo que o ótimo
de uma função objetivo nunca pode aumentar após um knockout.

A pergunta relevante para engenharia de cepas é outra, e é a mesma que o
OptKnock resolve: existe alguma deleção que torne a produção de CBDA
OBRIGATÓRIA quando a célula maximiza o próprio crescimento? A métrica
correspondente é o fluxo MÍNIMO de CBDA no ótimo de crescimento
(cbda_min > 0 indica acoplamento, isto é, estabilidade evolutiva da cepa).

Saídas: results/acoplamento_simples.csv, results/acoplamento_duplo.csv
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


def avaliar(modelo):
    """Retorna (mu, cbda_min, cbda_max) no ótimo de crescimento."""
    modelo.objective = BIO
    mu = modelo.slim_optimize()
    if mu is None or pd.isna(mu) or mu < TOL:
        return 0.0, 0.0, 0.0
    fva = flux_variability_analysis(
        modelo, reaction_list=[TARGET], fraction_of_optimum=0.999
    )
    return mu, float(fva.loc[TARGET, "minimum"]), float(fva.loc[TARGET, "maximum"])


def main():
    modelo = build()
    base = avaliar(modelo)
    print(f"selvagem: mu={base[0]:.4f} | CBDA no ótimo de crescimento: "
          f"min={base[1]:.4f} max={base[2]:.4f}")

    genes = [g.id for g in modelo.genes]
    linhas = []
    for i, gid in enumerate(genes):
        with modelo as mm:
            mm.genes.get_by_id(gid).knock_out()
            mu, cmin, cmax = avaliar(mm)
        linhas.append({"gene": gid, "mu": mu, "cbda_min": cmin, "cbda_max": cmax})
        if (i + 1) % 150 == 0:
            print(f"  ... {i + 1}/{len(genes)} genes")

    df = pd.DataFrame(linhas)
    df["mu_rel"] = df["mu"] / base[0]
    df["acoplado"] = df["cbda_min"] > TOL
    df = df.sort_values(["cbda_min", "mu"], ascending=False)
    df.to_csv(RES / "acoplamento_simples.csv", index=False)

    acop = df[df["acoplado"] & (df["mu_rel"] > 0.05)]
    print(f"\ndeleções simples com acoplamento (mu > 5% do selvagem): {len(acop)}")
    if len(acop):
        print(acop.head(20).to_string(index=False))

    # ------------------------------------------------------------------
    # Duplas: combina os genes que mais preservam crescimento, incluindo
    # alvos a priori citados no relatório parcial (ERG9, OLE1, glioxilato)
    # ------------------------------------------------------------------
    nomes_alvo = {"ERG9", "OLE1", "ICL1", "MLS1", "CIT2", "ERG20", "ADH1", "GPD1", "GPD2"}
    por_nome = {}
    for g in modelo.genes:
        rotulo = (g.name or "").upper()
        if rotulo in nomes_alvo:
            por_nome[rotulo] = g.id
    print(f"\nalvos a priori localizados no modelo: {por_nome}")

    viaveis = df[(df["mu_rel"] > 0.30)]["gene"].head(40).tolist()
    for gid in por_nome.values():
        if gid not in viaveis:
            viaveis.append(gid)

    linhas = []
    pares = list(combinations(viaveis, 2))
    print(f"avaliando {len(pares)} pares...")
    for a, b in pares:
        with modelo as mm:
            mm.genes.get_by_id(a).knock_out()
            mm.genes.get_by_id(b).knock_out()
            mu, cmin, cmax = avaliar(mm)
        linhas.append({"gene_a": a, "gene_b": b, "mu": mu,
                       "cbda_min": cmin, "cbda_max": cmax})

    dd = pd.DataFrame(linhas)
    dd["mu_rel"] = dd["mu"] / base[0]
    dd = dd.sort_values(["cbda_min", "mu"], ascending=False)
    dd.to_csv(RES / "acoplamento_duplo.csv", index=False)
    acop2 = dd[(dd["cbda_min"] > TOL) & (dd["mu_rel"] > 0.05)]
    print(f"\npares com acoplamento: {len(acop2)}")
    if len(acop2):
        print(acop2.head(20).to_string(index=False))


if __name__ == "__main__":
    main()
