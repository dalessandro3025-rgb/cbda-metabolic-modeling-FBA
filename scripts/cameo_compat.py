"""Shims para rodar cameo em Python moderno (importar ANTES de cameo)."""
import collections, collections.abc, pip
for _n in ("Mapping","MutableMapping","Iterable","Sequence","Callable"):
    if not hasattr(collections,_n): setattr(collections,_n,getattr(collections.abc,_n))
if not hasattr(pip,"get_installed_distributions"):
    pip.get_installed_distributions = lambda *a,**k: []
import sympy, sympy.core.singleton
if not hasattr(sympy, "singleton"):
    sympy.singleton = sympy.core.singleton

# cobra >= 0.26 removeu find_gene_knockout_reactions, ainda importada por cameo.
import cobra.manipulation.delete as _del
if not hasattr(_del, "find_gene_knockout_reactions"):
    def find_gene_knockout_reactions(model, gene_list, **kwargs):
        alvos = {g if isinstance(g, str) else g.id for g in gene_list}
        fora = []
        for r in model.reactions:
            if not r.genes or alvos.isdisjoint({g.id for g in r.genes}):
                continue
            if not r.gpr.eval(alvos):
                fora.append(r)
        return fora
    _del.find_gene_knockout_reactions = find_gene_knockout_reactions

# cameo usa pandas.core.common.in_ipython_frontend, movido nas versoes novas.
import pandas.core.common as _pdc
if not hasattr(_pdc, "in_ipython_frontend"):
    _pdc.in_ipython_frontend = lambda: False

# cobra >= 0.6 renomeou Solution.f para Solution.objective_value; cameo usa .f
from cobra.core.solution import Solution as _Sol
if not hasattr(_Sol, "f"):
    _Sol.f = property(lambda self: self.objective_value)
