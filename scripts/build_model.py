"""
Reconstrução do modelo iMM904 + via heteróloga do CBDA.

Centraliza em uma única função a inserção da via descrita no notebook
01_cbda_biosynthesis.ipynb, para que as análises subsequentes (FVA,
deleções gênicas, redox, OptKnock) partam exatamente do mesmo modelo.

Uso:
    from build_model import build
    modelo = build(hexanoate_exogenous=False)
"""

from pathlib import Path

import cobra
from cobra import Metabolite, Reaction

MODEL_PATH = Path(__file__).resolve().parent.parent / "models" / "iMM904.json"


def build(hexanoate_exogenous: bool = False, de_novo: bool = True, verbose: bool = False):
    """Carrega o iMM904 e insere a via heteróloga do CBDA.

    Parameters
    ----------
    hexanoate_exogenous : bool
        Adiciona a rota de ativação direta de hexanoato exógeno
        (EX_hexa_e -> hexanoil-CoA via acil-CoA ligase).
    de_novo : bool
        Mantém a rota de elongação de acetil-CoA (beta-oxidação reversa,
        módulos Clostridium + Euglena).
    """
    modelo = cobra.io.load_json_model(str(MODEL_PATH))

    # ------------------------------------------------------------------
    # Metabólitos heterólogos
    # ------------------------------------------------------------------
    novos = [
        Metabolite("butcoa_c", name="Butyryl-CoA", formula="C25H38N7O17P3S", charge=-4, compartment="c"),
        Metabolite("3hbutcoa_c", name="(S)-3-Hydroxybutyryl-CoA", formula="C25H38N7O18P3S", charge=-4, compartment="c"),
        Metabolite("crotcoa_c", name="Crotonyl-CoA", formula="C25H36N7O17P3S", charge=-4, compartment="c"),
        Metabolite("3ohexcoa_c", name="3-Oxohexanoyl-CoA", formula="C27H40N7O18P3S", charge=-4, compartment="c"),
        Metabolite("3hhexcoa_c", name="(S)-3-Hydroxyhexanoyl-CoA", formula="C27H42N7O18P3S", charge=-4, compartment="c"),
        Metabolite("t2hexcoa_c", name="trans-2-Hexenoyl-CoA", formula="C27H40N7O17P3S", charge=-4, compartment="c"),
        Metabolite("hexcoa_c", name="Hexanoyl-CoA", formula="C27H42N7O17P3S", charge=-4, compartment="c"),
        Metabolite("olivetolate_c", name="Olivetolate", formula="C12H15O4", charge=-1, compartment="c"),
        Metabolite("cbga_c", name="Cannabigerolate", formula="C22H31O4", charge=-1, compartment="c"),
        Metabolite("cbda_c", name="Cannabidiolate", formula="C22H29O4", charge=-1, compartment="c"),
    ]
    modelo.add_metabolites(novos)

    M = modelo.metabolites
    reacoes = []

    # ------------------------------------------------------------------
    # Rota A: elongação de novo (beta-oxidação reversa C4 -> C6)
    # ------------------------------------------------------------------
    if de_novo:
        r = Reaction("THL_c4", name="Thiolase (thlA) - C4 cycle")
        r.add_metabolites({M.accoa_c: -2, M.aacoa_c: 1, M.coa_c: 1})
        reacoes.append(r)

        r = Reaction("HBD_c4", name="3-hydroxyacyl-CoA dehydrogenase - C4")
        r.add_metabolites({M.aacoa_c: -1, M.nadh_c: -1, M.h_c: -1,
                           M.get_by_id("3hbutcoa_c"): 1, M.nad_c: 1})
        reacoes.append(r)

        r = Reaction("CRT_c4", name="Crotonase - C4")
        r.add_metabolites({M.get_by_id("3hbutcoa_c"): -1, M.crotcoa_c: 1, M.h2o_c: 1})
        reacoes.append(r)

        r = Reaction("TER_c4", name="trans-enoyl-CoA reductase - C4")
        r.add_metabolites({M.crotcoa_c: -1, M.nadh_c: -1, M.h_c: -1,
                           M.butcoa_c: 1, M.nad_c: 1})
        reacoes.append(r)

        r = Reaction("THL_c6", name="Thiolase (thlA) - C6 cycle")
        r.add_metabolites({M.butcoa_c: -1, M.accoa_c: -1,
                           M.get_by_id("3ohexcoa_c"): 1, M.coa_c: 1})
        reacoes.append(r)

        r = Reaction("HBD_c6", name="3-hydroxyacyl-CoA dehydrogenase - C6")
        r.add_metabolites({M.get_by_id("3ohexcoa_c"): -1, M.nadh_c: -1, M.h_c: -1,
                           M.get_by_id("3hhexcoa_c"): 1, M.nad_c: 1})
        reacoes.append(r)

        r = Reaction("CRT_c6", name="Crotonase - C6")
        r.add_metabolites({M.get_by_id("3hhexcoa_c"): -1, M.t2hexcoa_c: 1, M.h2o_c: 1})
        reacoes.append(r)

        r = Reaction("TER_c6", name="trans-enoyl-CoA reductase - C6")
        r.add_metabolites({M.t2hexcoa_c: -1, M.nadh_c: -1, M.h_c: -1,
                           M.hexcoa_c: 1, M.nad_c: 1})
        reacoes.append(r)

    # ------------------------------------------------------------------
    # Rota B: ativação de hexanoato exógeno
    # ------------------------------------------------------------------
    if hexanoate_exogenous:
        hexa_e = Metabolite("hexa_e", name="Hexanoate", formula="C6H11O2", charge=-1, compartment="e")
        hexa_c = Metabolite("hexa_c", name="Hexanoate", formula="C6H11O2", charge=-1, compartment="c")
        modelo.add_metabolites([hexa_e, hexa_c])

        r = Reaction("HEXAt", name="Hexanoate transport (diffusion)")
        r.add_metabolites({hexa_e: -1, hexa_c: 1})
        r.bounds = (-1000, 1000)
        reacoes.append(r)

        r = Reaction("HEXACL", name="Hexanoyl-CoA ligase (ACS-like)")
        r.add_metabolites({hexa_c: -1, M.coa_c: -1, M.atp_c: -1,
                           M.hexcoa_c: 1, M.amp_c: 1, M.ppi_c: 1})
        reacoes.append(r)

    # ------------------------------------------------------------------
    # Módulo Cannabis
    # ------------------------------------------------------------------
    r = Reaction("OLS_OAC", name="Tetraketide synthase / olivetolic acid cyclase")
    r.add_metabolites({M.hexcoa_c: -1, M.malcoa_c: -3, M.h_c: -2,
                       M.olivetolate_c: 1, M.co2_c: 3, M.coa_c: 4})
    reacoes.append(r)

    r = Reaction("CBGAS", name="Prenyltransferase CsPT4")
    r.add_metabolites({M.olivetolate_c: -1, M.grdp_c: -1, M.cbga_c: 1, M.ppi_c: 1})
    reacoes.append(r)

    r = Reaction("CBDAS", name="CBDA synthase (FAD-dependent oxidase)")
    r.add_metabolites({M.cbga_c: -1, M.o2_c: -1, M.cbda_c: 1, M.h2o2_c: 1})
    reacoes.append(r)

    modelo.add_reactions(reacoes)

    modelo.add_boundary(M.cbda_c, type="demand", reaction_id="DM_cbda_c")
    if hexanoate_exogenous:
        modelo.add_boundary(modelo.metabolites.hexa_e, type="exchange",
                            reaction_id="EX_hexa_e", lb=-1000, ub=1000)

    if verbose:
        for rid in [rx.id for rx in reacoes]:
            desvio = modelo.reactions.get_by_id(rid).check_mass_balance()
            if desvio:
                print(f"[balanço] {rid}: {desvio}")

    return modelo


if __name__ == "__main__":
    m = build(verbose=True)
    m.objective = "DM_cbda_c"
    s = m.optimize()
    print(f"\nreações: {len(m.reactions)} | metabólitos: {len(m.metabolites)}")
    print(f"status: {s.status} | CBDA máximo: {s.objective_value:.4f} mmol/gDCW/h")

    m2 = build(hexanoate_exogenous=True)
    m2.objective = "DM_cbda_c"
    print(f"com hexanoato exógeno: {m2.optimize().objective_value:.4f} mmol/gDCW/h")
