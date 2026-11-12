"""
Etapa 4: o que de fato limita o CBDA — precursores e balanço redox.

Estratégia: em vez de deleções (incapazes de elevar o ótimo), quantifica-se
a sensibilidade do fluxo de CBDA a (i) preços-sombra dos metabólitos-chave,
(ii) suplementação artificial de cada precursor, e (iii) intervenções de
balanço redox (transidrogenase heteróloga, via das pentoses fosfato).
"""

import sys
from pathlib import Path

import pandas as pd
from cobra import Reaction

sys.path.insert(0, str(Path(__file__).resolve().parent))
from build_model import build

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

TARGET = "DM_cbda_c"
CHAVE = ["accoa_c", "malcoa_c", "hexcoa_c", "grdp_c", "ipdp_c", "nadph_c",
         "nadh_c", "atp_c", "o2_c", "olivetolate_c", "cbga_c"]


def base_model():
    m = build()
    m.objective = TARGET
    return m


def precos_sombra():
    m = base_model()
    sol = m.optimize()
    linhas = []
    for mid in CHAVE:
        try:
            met = m.metabolites.get_by_id(mid)
        except KeyError:
            continue
        linhas.append({"metabolito": mid, "nome": met.name,
                       "preco_sombra": float(sol.shadow_prices.get(mid, float("nan")))})
    return pd.DataFrame(linhas)


def suplementacao():
    """Adiciona um sink de entrada (1 mmol/gDCW/h) para cada precursor e mede o ganho."""
    ref = base_model().slim_optimize()
    linhas = []
    for mid in CHAVE:
        m = base_model()
        try:
            met = m.metabolites.get_by_id(mid)
        except KeyError:
            continue
        r = Reaction(f"SUP_{mid}")
        r.add_metabolites({met: 1})
        r.bounds = (0, 1)
        m.add_reactions([r])
        v = m.slim_optimize()
        linhas.append({"precursor": mid, "cbda": v, "ganho_pct": 100 * (v - ref) / ref})
    return ref, pd.DataFrame(linhas).sort_values("ganho_pct", ascending=False)


def redox():
    ref = base_model().slim_optimize()
    linhas = [{"cenario": "referência", "cbda": ref, "delta_pct": 0.0}]

    # (a) transidrogenase heteróloga NADH -> NADPH
    m = base_model()
    r = Reaction("THD_het", name="Transidrogenase heteróloga (NADH -> NADPH)")
    r.add_metabolites({m.metabolites.nadh_c: -1, m.metabolites.nadp_c: -1,
                       m.metabolites.nad_c: 1, m.metabolites.nadph_c: 1})
    r.bounds = (0, 1000)
    m.add_reactions([r])
    v = m.slim_optimize()
    linhas.append({"cenario": "transidrogenase NADH->NADPH", "cbda": v,
                   "delta_pct": 100 * (v - ref) / ref})

    # (b) bloqueio da via das pentoses fosfato (G6PDH2r / ZWF1)
    m = base_model()
    for rid in ["G6PDH2r", "G6PDH2"]:
        if rid in [x.id for x in m.reactions]:
            m.reactions.get_by_id(rid).bounds = (0, 0)
    v = m.slim_optimize()
    linhas.append({"cenario": "PPP oxidativa bloqueada", "cbda": v,
                   "delta_pct": 100 * (v - ref) / ref})

    # (c) ERG9 (esqualeno sintase) bloqueada — dreno de FPP
    m = base_model()
    g = next((x for x in m.genes if (x.name or "").upper() == "ERG9"), None)
    if g:
        with m:
            g.knock_out()
            v = m.slim_optimize()
        linhas.append({"cenario": "ERG9 (esqualeno sintase) deletada", "cbda": v,
                       "delta_pct": 100 * (v - ref) / ref})

    # (d) ATP de manutenção
    for atpm in ["ATPM"]:
        if atpm in [x.id for x in base_model().reactions]:
            for fator, rotulo in [(0.0, "ATPM = 0"), (2.0, "ATPM dobrado")]:
                m = base_model()
                rx = m.reactions.get_by_id(atpm)
                alvo = rx.lower_bound * fator
                rx.bounds = (min(alvo, rx.upper_bound), max(alvo, rx.upper_bound))
                rx.bounds = (alvo, max(alvo, rx.upper_bound))
                v = m.slim_optimize()
                linhas.append({"cenario": rotulo, "cbda": v,
                               "delta_pct": 100 * (v - ref) / ref})

    # (e) O2 irrestrito
    m = base_model()
    m.reactions.EX_o2_e.lower_bound = -1000
    v = m.slim_optimize()
    linhas.append({"cenario": "O2 irrestrito", "cbda": v, "delta_pct": 100 * (v - ref) / ref})

    return pd.DataFrame(linhas)


if __name__ == "__main__":
    ps = precos_sombra()
    ps.to_csv(RES / "precos_sombra.csv", index=False)
    print("=== preços-sombra no ótimo de CBDA ===")
    print(ps.to_string(index=False))

    ref, sup = suplementacao()
    sup.to_csv(RES / "suplementacao_precursores.csv", index=False)
    print(f"\n=== suplementação de 1 mmol/gDCW/h (referência = {ref:.4f}) ===")
    print(sup.to_string(index=False))

    rx = redox()
    rx.to_csv(RES / "cenarios_redox.csv", index=False)
    print("\n=== cenários de balanço redox e energia ===")
    print(rx.to_string(index=False))
