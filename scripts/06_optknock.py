"""
Etapa 6: OptKnock formal (cameo), para registro do resultado no relatório.

Ambiente: cameo 0.11.2 não roda com cobra >= 0.27 nem numpy >= 2. Este
script exige o ambiente isolado descrito em docs/ambiente_cameo.md
(Python 3.11, numpy < 2, cobra 0.26.3, cameo) e importa cameo_compat
ANTES de cameo para aplicar os shims de compatibilidade necessários.

Execução:
    ../venv_cameo/bin/python scripts/06_optknock.py
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cameo_compat  # noqa: F401  (aplica shims; deve preceder cameo)

from cameo.strain_design.deterministic.linear_programming import OptKnock  # noqa: E402

from build_model import build  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)

BIO = "BIOMASS_SC5_notrace"
TARGET = "DM_cbda_c"

SUBSISTEMAS = ("glycolysis", "pentose", "citric", "pyruvate", "oxidative",
               "anaplerotic", "sterol", "steroid", "fatty acid", "glyoxylate",
               "alternative carbon")


def espaco_de_busca(modelo):
    """Reações deletáveis: têm GPR, não são de fronteira, no metabolismo central."""
    sel = []
    for r in modelo.reactions:
        if r.id in (BIO, TARGET) or r.boundary or not r.genes:
            continue
        sub = (r.subsystem or "").lower()
        if any(s in sub for s in SUBSISTEMAS):
            sel.append(r.id)
    return sel


def main():
    modelo = build()
    modelo.objective = BIO

    alvos = espaco_de_busca(modelo)
    nao_alvos = [r.id for r in modelo.reactions if r.id not in alvos]
    print(f"reações no modelo: {len(modelo.reactions)}")
    print(f"espaço de busca do OptKnock: {len(alvos)} reações deletáveis", flush=True)

    ok = OptKnock(modelo, exclude_reactions=nao_alvos, remove_blocked=False,
                  fraction_of_optimum=0.1)

    resumo = []
    for max_ko in (1, 2, 3):
        print(f"\n--- OptKnock, ate {max_ko} delecao(oes) ---", flush=True)
        try:
            res = ok.run(max_knockouts=max_ko, target=TARGET, biomass=BIO,
                         max_results=5)
            df = res.data_frame
            n = len(df)
            print(f"solucoes: {n}")
            if n:
                print(df.to_string())
                df.to_csv(RES / f"optknock_{max_ko}ko.csv", index=False)
            resumo.append((max_ko, n))
        except Exception as exc:
            print(f"falha: {type(exc).__name__}: {exc}")
            resumo.append((max_ko, None))

    print("\n=== resumo ===")
    for k, n in resumo:
        print(f"  ate {k} delecao(oes): {n if n is not None else 'erro'} solucao(oes)")


if __name__ == "__main__":
    main()
