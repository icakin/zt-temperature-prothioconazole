#!/usr/bin/env python3
# build_gecko_full.py — FULL GECKO expansion of ztGEM v0.3:
# isozymes (OR arms) become separate reaction arms; complexes (AND) consume
# every subunit's prot_<gene>; reversible arms split; draw reactions per gene
# with real MWs. kcat remains the uniform placeholder (DLTKcat pending).
import cobra, pandas as pd, numpy as np, logging, re, ast, itertools
from Bio import SeqIO
logging.getLogger("cobra").setLevel(logging.ERROR)

KCAT_H = 25.0*3600.0
POOL = 0.5*0.4*0.45
MAX_ARMS = 12

aa_w = {'A':71.08,'R':156.19,'N':114.10,'D':115.09,'C':103.14,'E':129.12,'Q':128.13,
        'G':57.05,'H':137.14,'I':113.16,'L':113.16,'K':128.17,'M':131.19,'F':147.18,
        'P':97.12,'S':87.08,'T':101.10,'W':186.21,'Y':163.18,'V':99.13}
mw = {}
for rec in SeqIO.parse("Zymoseptoria_tritici.MG2.pep.all.fa","fasta"):
    m = re.search(r"gene:(\S+)", rec.description)
    if not m: continue
    g = m.group(1)
    w = sum(aa_w.get(a,110.0) for a in str(rec.seq).rstrip("*")) + 18.02
    mw[g] = max(mw.get(g,0), w/1000.0)

M = cobra.io.read_sbml_model("models/gem/ztGEM_v03.xml")
M.reactions.r_4046.bounds = (0.0, 1000.0)

def gpr_arms(rule):
    """Return list of complexes (each = list of Zt genes) from a GPR string."""
    if not rule: return []
    try:
        tree = ast.parse(rule.replace(" and ", " & ").replace(" or ", " | ")
                              .replace("and", "&").replace("or", "|"), mode="eval")
    except SyntaxError:
        return []
    def walk(n):
        if isinstance(n, ast.Name): return [[n.id]]
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.BitOr):
            return walk(n.left) + walk(n.right)
        if isinstance(n, ast.BinOp) and isinstance(n.op, ast.BitAnd):
            L, R = walk(n.left), walk(n.right)
            out = [l + r for l, r in itertools.product(L, R)]
            return out[:MAX_ARMS]
        if isinstance(n, ast.Expr): return walk(n.value)
        return []
    arms = walk(tree.body)
    # keep Zt subunits only; drop empty
    arms = [[g for g in a if g.startswith("Mycgr3G")] for a in arms]
    arms = [a for a in arms if a]
    # dedupe
    seen = set(); out = []
    for a in arms:
        k = tuple(sorted(set(a)))
        if k not in seen: seen.add(k); out.append(sorted(set(a)))
    return out[:MAX_ARMS]

cyto = "c"
prot_mets = {}
def pm(g):
    if g not in prot_mets:
        prot_mets[g] = cobra.Metabolite(f"prot_{g}", name=f"enzyme {g}", compartment=cyto)
    return prot_mets[g]

new, to_remove = [], []
n_arms_total = 0
skip = {"r_2111","r_4041","r_4046","r_4048","r_4049","r_4050"}
for rxn in list(M.reactions):
    if rxn.boundary or rxn.id in skip: continue
    arms = gpr_arms(rxn.gene_reaction_rule)
    if not arms: continue
    directions = []
    if rxn.upper_bound > 0: directions.append(("", 1.0, rxn.upper_bound))
    if rxn.lower_bound < 0: directions.append(("_REV", -1.0, abs(rxn.lower_bound)))
    for ai, complex_genes in enumerate(arms, 1):
        for tag, sign, ub in directions:
            a = cobra.Reaction(f"{rxn.id}_No{ai}{tag}",
                               name=f"{rxn.name or rxn.id} (isozyme {ai}{' rev' if tag else ''})")
            a.add_metabolites({met: sign*c for met, c in rxn.metabolites.items()})
            a.bounds = (0, ub)
            a.gene_reaction_rule = " and ".join(complex_genes)
            a.add_metabolites({pm(g): -1.0/KCAT_H for g in complex_genes})
            new.append(a); n_arms_total += 1
    to_remove.append(rxn)
M.add_reactions(new)
M.remove_reactions(to_remove)
print(f"expanded {len(to_remove)} reactions into {n_arms_total} enzyme arms; {len(prot_mets)} enzymes")

pool = cobra.Metabolite("prot_pool", name="protein pool", compartment=cyto)
draws = []
for g, met in prot_mets.items():
    d = cobra.Reaction(f"draw_prot_{g}")
    d.add_metabolites({pool: -mw.get(g, 0.045), met: 1.0}); d.bounds = (0, 1000)
    draws.append(d)
pex = cobra.Reaction("prot_pool_exchange"); pex.add_metabolites({pool: 1.0})
pex.bounds = (0, POOL)
M.add_reactions(draws + [pex])

for g in list(M.groups):
    if g.name: g.name = re.sub(r"[\s,/]+", "_", g.name)

med = M.medium
for rid in ["r_1714","r_4502","r_4504"]:
    if rid in med: med[rid] = 0
med["r_2058"] = 1000
M.medium = med
mu = M.slim_optimize()
sol = M.optimize()
print("ecZtGEM-full growth:", round(mu,4), " pool use:", round(sol.fluxes["prot_pool_exchange"],4), "/", POOL)
M.id = "ecZtGEM_full"
M.name = "Z. tritici IPO323 enzyme-constrained GEM (full GECKO expansion; placeholder kcats)"
cobra.io.write_sbml_model(M, "models/gem/ecZtGEM_full.xml")
print("saved ecZtGEM_full.xml:", len(M.reactions), "rxns,", len(M.metabolites), "mets")
