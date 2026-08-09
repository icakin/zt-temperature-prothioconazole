#!/usr/bin/env python3
# build_ztGEM_v03.py — ztGEM v0.3: diamond-BBH orthology, Zt GPRs, AOX, gap-fill.
import cobra, pandas as pd, json, logging, re
logging.getLogger("cobra").setLevel(logging.ERROR)

Y = cobra.io.read_sbml_model("models/gem/scaffolds/yeast-GEM_v9.1.0.xml")
bbh = pd.read_csv("models/gem/bbh_zt_sc.csv")
sc2zt = {}
for _, r in bbh.iterrows():
    sc2zt.setdefault(r["sc_gene"], []).append(r["zt_gene"])

# one-way rescue (yeast->Zt best hit, decent score) for enzymes missed by BBH
cols = ["q","s","pid","len","mm","go","qs","qe","ss","se","ev","bit"]
ow = pd.read_csv("models/gem/sc2zt.tsv", sep="\t", names=cols).sort_values("bit", ascending=False).drop_duplicates("q")
ow = ow[(ow["bit"] >= 100)]
ow_map = {r["q"]: r["s"].replace("Mycgr3P","Mycgr3G") for _, r in ow.iterrows()}

anno = pd.read_csv("data/reference/gene_annotation.csv")
anno["EC"] = anno["EC"].fillna("")
zt_ecs = set()
for e in anno["EC"]: zt_ecs.update(x for x in str(e).split(";") if x and x != "nan")

def rxn_ecs(rxn):
    ec = rxn.annotation.get("ec-code", [])
    return set([ec] if isinstance(ec, str) else ec)

keep = {}
drop = []
for rxn in Y.reactions:
    genes = list(rxn.genes)
    if not genes:
        keep[rxn.id] = "no-gene"; continue
    if any(g.id in sc2zt for g in genes):
        keep[rxn.id] = "BBH"; continue
    if any(g.id in ow_map for g in genes):
        keep[rxn.id] = "one-way"; continue
    if rxn_ecs(rxn) & zt_ecs:
        keep[rxn.id] = "EC"; continue
    comps = {m.compartment for m in rxn.metabolites}
    if len(comps) > 1 and not rxn_ecs(rxn):
        keep[rxn.id] = "transport-flagged"; continue
    drop.append(rxn.id)
from collections import Counter
print(Counter(keep.values()), "drop:", len(drop))

M = Y.copy()
M.remove_reactions([M.reactions.get_by_id(r) for r in drop], remove_orphans=True)

## rewrite GPRs to Zt genes where orthologs exist
def zt_of(scg):
    if scg in sc2zt: return sc2zt[scg]
    if scg in ow_map: return [ow_map[scg]]
    return None
rewritten = 0
for rxn in M.reactions:
    if not rxn.genes: continue
    gpr = rxn.gene_reaction_rule
    parts = sorted({g.id for g in rxn.genes}, key=len, reverse=True)
    new = gpr
    ok = True
    for g in parts:
        z = zt_of(g)
        if z is None:
            ok = False; continue
        new = re.sub(rf"\b{g}\b", "( " + " or ".join(z) + " )" if len(z) > 1 else z[0], new)
    rxn.gene_reaction_rule = new
    if ok: rewritten += 1
print("GPRs fully rewritten to Zt genes:", rewritten)

## AOX
uqh2 = M.metabolites.get_by_id("s_1535"); uq = M.metabolites.get_by_id("s_1537")
o2m = M.metabolites.get_by_id("s_1278"); h2om = M.metabolites.get_by_id("s_0807")
aox = cobra.Reaction("r_AOX"); aox.name = "alternative oxidase (non-phosphorylating)"
aox.add_metabolites({uqh2: -1, o2m: -0.5, uq: 1, h2om: 1}); aox.bounds = (0, 1000)
aox.gene_reaction_rule = "Mycgr3G72918"
M.add_reactions([aox])

## medium: sucrose
med = M.medium
for rid in ["r_1714", "r_4502", "r_4504"]:
    if rid in med: med[rid] = 0
med["r_2058"] = 10
M.medium = med

## gap-fill greedily from dropped set
F = M.copy()
F.add_reactions([Y.reactions.get_by_id(i).copy() for i in drop if i not in [r.id for r in F.reactions]])
med = F.medium
for rid in ["r_1714","r_4502","r_4504"]:
    if rid in med: med[rid] = 0
med["r_2058"] = 10
F.medium = med
base = F.slim_optimize(); print("full growth:", round(base, 4))
kept_back = []
for rid in sorted(drop):
    r = F.reactions.get_by_id(rid); lb, ub = r.bounds
    r.bounds = (0, 0)
    v = F.slim_optimize()
    if v is None or v != v or v < 0.5 * base:
        r.bounds = (lb, ub); kept_back.append(rid)
F.remove_reactions([F.reactions.get_by_id(i) for i in sorted(drop) if i not in kept_back],
                   remove_orphans=True)
print("gap-fill add-backs:", [(i, Y.reactions.get_by_id(i).name[:50]) for i in kept_back])
print("v0.3:", len(F.reactions), "rxns; growth on sucrose:", round(F.slim_optimize(), 4))

F.id = "ztGEM_v03"
F.name = "Z. tritici IPO323 GEM v0.3 (Yeast9 template, diamond BBH orthology, AOX added)"
cobra.io.write_sbml_model(F, "models/gem/ztGEM_v03.xml")
json.dump({"support": keep, "gapfill": kept_back}, open("models/gem/ztGEM_v03_provenance.json", "w"))
print("saved ztGEM_v03.xml")

## key checks
erg = [m for m in F.metabolites if (m.name or "").lower() == "ergosterol"]
cyp = F.reactions.get_by_id("r_0317")
print("CYP51 GPR now:", cyp.gene_reaction_rule)
print("AOX GPR:", F.reactions.r_AOX.gene_reaction_rule)
