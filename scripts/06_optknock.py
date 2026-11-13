"""
OptKnock formal (cameo), resiliente a parede de tempo.

Executa a busca em estagios crescentes (1, 2, 3 delecoes) e grava o
resultado de cada estagio assim que ele termina, em
results/optknock_status.json. Se o job cair durante o estagio de 3
delecoes, os estagios de 1 e 2 ja estarao registrados.

Execucao:  python -u scripts/06_optknock.py [max_results]
"""

import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import cameo_compat  # noqa: F401  (shims; precede cameo)

from cameo.strain_design.deterministic.linear_programming import OptKnock  # noqa: E402

from build_model import build  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
RES.mkdir(exist_ok=True)
STATUS = RES / "optknock_status.json"

BIO = "BIOMASS_SC5_notrace"
TARGET = "DM_cbda_c"

SUBSISTEMAS = ("glycolysis", "pentose", "citric", "pyruvate", "oxidative",
               "anaplerotic", "sterol", "steroid", "fatty acid", "glyoxylate",
               "alternative carbon")


def agora():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def gravar(estado):
    STATUS.write_text(json.dumps(estado, indent=2, ensure_ascii=False))


def espaco_de_busca(modelo):
    sel = []
    for r in modelo.reactions:
        if r.id in (BIO, TARGET) or r.boundary or not r.genes:
            continue
        if any(s in (r.subsystem or "").lower() for s in SUBSISTEMAS):
            sel.append(r.id)
    return sel


def main():
    max_results = int(sys.argv[1]) if len(sys.argv) > 1 else 3

    estado = {"inicio": agora(), "max_results": max_results, "estagios": []}
    gravar(estado)

    modelo = build()
    modelo.objective = BIO

    alvos = espaco_de_busca(modelo)
    nao_alvos = [r.id for r in modelo.reactions if r.id not in alvos]
    estado["espaco_de_busca"] = len(alvos)
    estado["reacoes_no_modelo"] = len(modelo.reactions)
    gravar(estado)

    print(f"[{agora()}] reacoes no modelo: {len(modelo.reactions)}", flush=True)
    print(f"[{agora()}] espaco de busca: {len(alvos)} reacoes deletaveis", flush=True)

    t0 = time.time()
    ok = OptKnock(modelo, exclude_reactions=nao_alvos, remove_blocked=False,
                  fraction_of_optimum=0.1)
    estado["montagem_do_milp_s"] = round(time.time() - t0, 1)
    gravar(estado)
    print(f"[{agora()}] MILP montado em {estado['montagem_do_milp_s']}s", flush=True)

    for max_ko in (1, 2, 3):
        print(f"\n[{agora()}] ===== ate {max_ko} delecao(oes) =====", flush=True)
        estagio = {"max_knockouts": max_ko, "inicio": agora()}
        estado["estagios"].append(estagio)
        gravar(estado)

        t0 = time.time()
        try:
            res = ok.run(max_knockouts=max_ko, target=TARGET, biomass=BIO,
                         max_results=max_results)
            df = res.data_frame
            estagio["solucoes"] = len(df)
            estagio["status"] = "concluido"
            if len(df):
                destino = RES / f"optknock_{max_ko}ko.csv"
                df.to_csv(destino, index=False)
                estagio["arquivo"] = destino.name
                print(df.to_string(), flush=True)
            else:
                print("nenhum conjunto de delecoes acopla crescimento a CBDA.",
                      flush=True)
        except Exception as exc:
            estagio["status"] = f"{type(exc).__name__}: {exc}"
            print(f"falha: {estagio['status']}", flush=True)

        estagio["duracao_s"] = round(time.time() - t0, 1)
        estagio["fim"] = agora()
        gravar(estado)
        print(f"[{agora()}] estagio {max_ko} encerrado em "
              f"{estagio['duracao_s']}s: {estagio['status']}", flush=True)

    estado["fim"] = agora()
    estado["completo"] = True
    gravar(estado)

    print(f"\n[{agora()}] ===== resumo =====", flush=True)
    for e in estado["estagios"]:
        print(f"  ate {e['max_knockouts']} KO: {e.get('status')} "
              f"({e.get('solucoes', '-')} solucoes, {e['duracao_s']}s)", flush=True)


if __name__ == "__main__":
    main()
