"""
Etapa 2 do projeto: envelope de produção, FVA e deleções gênicas.

Motivação: a linha de base de 0,7885 mmol/gDCW/h do relatório parcial
corresponde ao ótimo de CBDA com crescimento NULO. Para que as análises
de engenharia tenham significado fisiológico, todas as simulações aqui
impõem um crescimento mínimo (fração de mu_max).

Saídas: results/*.csv e figures/*.png
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cobra.flux_analysis import (
    flux_variability_analysis,
    single_gene_deletion,
    double_gene_deletion,
    production_envelope,
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_model import build

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
FIG = ROOT / "figures"
RES.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)

BIO = "BIOMASS_SC5_notrace"
TARGET = "DM_cbda_c"
GROWTH_FRACTION = 0.30  # crescimento mínimo imposto, fração de mu_max


def main():
    modelo = build()

    modelo.objective = BIO
    mu_max = modelo.slim_optimize()
    modelo.objective = TARGET
    cbda_max = modelo.slim_optimize()
    print(f"mu_max = {mu_max:.4f} /h | CBDA_max (crescimento nulo) = {cbda_max:.4f} mmol/gDCW/h")

    # ------------------------------------------------------------------
    # 1. Envelope de produção (trade-off crescimento x CBDA)
    # ------------------------------------------------------------------
    env = production_envelope(modelo, [BIO], objective=TARGET, points=25)
    env.to_csv(RES / "production_envelope.csv", index=False)

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.plot(env[BIO], env["flux_maximum"], marker="o", ms=4, label="CBDA máximo")
    ax.plot(env[BIO], env["flux_minimum"], ls="--", label="CBDA mínimo")
    ax.set_xlabel("Taxa específica de crescimento (h$^{-1}$)")
    ax.set_ylabel("Fluxo de CBDA (mmol gDCW$^{-1}$ h$^{-1}$)")
    ax.set_title("Envelope de produção: trade-off crescimento x CBDA")
    ax.legend()
    ax.grid(alpha=.3)
    fig.tight_layout()
    fig.savefig(FIG / "production_envelope.png", dpi=300)
    plt.close(fig)
    print(f"envelope salvo ({len(env)} pontos)")

    # ------------------------------------------------------------------
    # Cenário de referência: crescimento mínimo imposto
    # ------------------------------------------------------------------
    mu_min = GROWTH_FRACTION * mu_max
    modelo.reactions.get_by_id(BIO).lower_bound = mu_min
    modelo.objective = TARGET
    cbda_ref = modelo.slim_optimize()
    print(f"cenário de referência: mu >= {mu_min:.4f} /h -> CBDA = {cbda_ref:.4f} mmol/gDCW/h")

    # ------------------------------------------------------------------
    # 2. FVA
    # ------------------------------------------------------------------
    fva = flux_variability_analysis(modelo, fraction_of_optimum=0.95)
    fva["amplitude"] = fva["maximum"] - fva["minimum"]
    fva["essencial"] = (fva["amplitude"] < 1e-6) & (fva["minimum"].abs() > 1e-6)
    fva.sort_values("amplitude", ascending=False).to_csv(RES / "fva.csv")
    print(f"FVA: {int(fva['essencial'].sum())} reações de fluxo fixo não nulo "
          f"(essenciais p/ o ótimo); {int((fva['amplitude'] > 1).sum())} com amplitude > 1")

    # ------------------------------------------------------------------
    # 3. Deleção gênica simples
    # ------------------------------------------------------------------
    sgd = single_gene_deletion(modelo, processes=1)
    sgd["gene"] = [list(ids)[0] for ids in sgd.index.get_level_values(0)] \
        if isinstance(sgd.index, pd.MultiIndex) else [list(i)[0] for i in sgd["ids"]]
    sgd = sgd.reset_index(drop=True)[["gene", "growth", "status"]]
    sgd = sgd.rename(columns={"growth": "cbda"})
    sgd["delta_pct"] = 100 * (sgd["cbda"] - cbda_ref) / cbda_ref
    sgd = sgd.sort_values("cbda", ascending=False)
    sgd.to_csv(RES / "single_gene_deletion.csv", index=False)

    melhoram = sgd[sgd["delta_pct"] > 1e-6]
    letais = sgd[sgd["cbda"].fillna(0) < 1e-6]
    print(f"deleção simples: {len(sgd)} genes | {len(melhoram)} aumentam CBDA | "
          f"{len(letais)} inviabilizam o cenário")
    if len(melhoram):
        print(melhoram.head(15).to_string(index=False))

    # ------------------------------------------------------------------
    # 4. Deleção dupla sobre os melhores candidatos
    # ------------------------------------------------------------------
    cand = list(melhoram["gene"].head(25))
    if len(cand) >= 2:
        dgd = double_gene_deletion(modelo, gene_list1=cand, processes=1)
        dgd = dgd.reset_index(drop=True)
        dgd["genes"] = dgd["ids"].apply(lambda s: " + ".join(sorted(s)))
        dgd = dgd.rename(columns={"growth": "cbda"})[["genes", "cbda", "status"]]
        dgd["delta_pct"] = 100 * (dgd["cbda"] - cbda_ref) / cbda_ref
        dgd = dgd.sort_values("cbda", ascending=False)
        dgd.to_csv(RES / "double_gene_deletion.csv", index=False)
        print("\ndeleção dupla (top 10):")
        print(dgd.head(10).to_string(index=False))
    else:
        print("candidatos insuficientes para deleção dupla")


if __name__ == "__main__":
    main()
