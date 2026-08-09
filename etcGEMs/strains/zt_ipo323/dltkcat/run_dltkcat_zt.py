#!/usr/bin/env python3
"""run_dltkcat_zt.py — run DLTKcat temperature-dependent kcat predictions for
ztGEM (zt_ipo323 strain) without DLTKcat's web lookups.

Replaces convert_input (UniProt/PubChem web calls) with local sources:
  - sequences: Zymoseptoria_tritici.MG2.pep.all.fa (gene: tags)
  - SMILES:    dltkcat_subs.csv (ModelSEED compounds.tsv via BiGG id + training-name match)
Replicates gen_features + predict.py exactly, with one deliberate fix: unknown
dict items (atoms/fingerprints/words unseen in training) map to the single
spare embedding slot (index = len(dict)) instead of growing the dict, which
would overflow the n+1-sized embeddings.

Temperature normalisation uses the training-set min/max from
data/processed_data.csv: Temp_K in [273.15, 373.15].

Output: dltkcat_output.csv with rxn_id, enz, sub, Temp_C, pred_log10kcat
(= strains/zt_ipo323/dltkcat/output.csv for `etcgem dltkcat parse`).
"""
import os, sys, pickle, math
import numpy as np
import pandas as pd
import torch
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit import RDLogger
RDLogger.DisableLog("rdApp.*")

HERE = os.path.dirname(os.path.abspath(__file__))
DLT = os.path.join(HERE, "DLTKcat")
sys.path.insert(0, os.path.join(DLT, "code"))
from DLTKcat import DLTKcat  # noqa: E402

INPUT = os.path.join(HERE, "dltkcat_input_enriched.csv")  # smiles pre-attached
FASTA = os.path.join(HERE, "Zymoseptoria_tritici.MG2.pep.all.fa")
OUT = os.path.join(HERE, "dltkcat_output.csv")
MODEL = os.path.join(DLT, "data", "performances",
                     "model_latentdim=40_outlayer=4_rmsetest=0.8854_rmsedev=0.908.pth")

TK_MIN, TK_MAX = 273.15, 373.15
INV_MIN, INV_MAX = 1.0 / TK_MAX, 1.0 / TK_MIN
RADIUS, NGRAM = 2, 3

def load_pkl(p):
    with open(p, "rb") as f:
        return pickle.load(f)

dict_dir = os.path.join(DLT, "data", "dict")
atom_dict = load_pkl(os.path.join(dict_dir, "atom_dict.pkl"))
bond_dict = load_pkl(os.path.join(dict_dir, "bond_dict.pkl"))
fp_dict = load_pkl(os.path.join(dict_dir, "fingerprint_dict.pkl"))
edge_dict = load_pkl(os.path.join(dict_dir, "edge_dict.pkl"))
word_dict = load_pkl(os.path.join(dict_dir, "word_dict.pkl"))
SPARE = {id(atom_dict): len(atom_dict), id(bond_dict): len(bond_dict),
         id(fp_dict): len(fp_dict), id(edge_dict): len(edge_dict),
         id(word_dict): len(word_dict)}

def look(item, d):
    """frozen check_dict: unseen -> spare slot len(d) (embeddings are n+1)."""
    return d.get(item, SPARE[id(d)])

def compound_features(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    mol = Chem.AddHs(mol)
    if mol.GetNumAtoms() < 2:
        return None
    atoms = [a.GetSymbol() for a in mol.GetAtoms()]
    for a in mol.GetAromaticAtoms():
        atoms[a.GetIdx()] = (atoms[a.GetIdx()], "aromatic")
    atoms = np.array([look(a, atom_dict) for a in atoms])
    ijb = defaultdict(list)
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        bd = look(str(b.GetBondType()), bond_dict)
        ijb[i].append((j, bd)); ijb[j].append((i, bd))
    iso = set(range(mol.GetNumAtoms())) - set(ijb.keys())
    for a in iso:
        ijb[a].append((a, look("nan", bond_dict)))
    # radius-2 fingerprints
    nodes, ijedge = atoms, ijb
    fingerprints = None
    for _ in range(RADIUS):
        fingerprints = []
        for i, j_edge in ijedge.items():
            neighbors = [(nodes[j], e) for j, e in j_edge]
            fingerprints.append(look((nodes[i], tuple(sorted(neighbors))), fp_dict))
        nodes = fingerprints
        nxt = defaultdict(list)
        for i, j_edge in ijedge.items():
            for j, e in j_edge:
                both = tuple(sorted((nodes[i], nodes[j])))
                nxt[i].append((j, look((both, e), edge_dict)))
        ijedge = nxt
    fp_arr = np.zeros((0,), dtype=np.int8)
    bi = {}
    fp = AllChem.GetMorganFingerprintAsBitVect(mol, RADIUS, nBits=1024,
                                               useChirality=True, bitInfo=bi)
    DataStructs.ConvertToNumpyArray(fp, fp_arr)
    adj = np.array(Chem.GetAdjacencyMatrix(mol)) + np.eye(mol.GetNumAtoms(), dtype=int)
    return np.array(fingerprints), adj, fp_arr

def protein_words(seq):
    s = ">" + seq + "<"
    return np.array([look(s[i:i + NGRAM], word_dict) for i in range(len(s) - NGRAM + 1)])

# ---------- assemble table ----------
inp = pd.read_csv(INPUT)

seqs = {}
gid, cur = None, []
for line in open(FASTA):
    if line.startswith(">"):
        if gid and cur:
            s = "".join(cur).rstrip("*")
            if len(s) > len(seqs.get(gid, "")):
                seqs[gid] = s
        import re
        m = re.search(r"gene:(\S+)", line)
        gid = m.group(1) if m else None
        cur = []
    else:
        cur.append(line.strip())
if gid and cur:
    s = "".join(cur).rstrip("*")
    if len(s) > len(seqs.get(gid, "")):
        seqs[gid] = s
inp["seq"] = inp["enz"].map(seqs)

n0 = len(inp)
inp = inp.dropna(subset=["smiles", "seq"]).reset_index(drop=True)
print(f"rows {n0} -> {len(inp)} after smiles/seq attach "
      f"({inp['rxn_id'].nunique()} reactions)", flush=True)

# unique compound / protein features
comp_cache = {}
bad_smiles = set()
for sm in inp["smiles"].unique():
    f = compound_features(sm)
    if f is None:
        bad_smiles.add(sm)
    else:
        comp_cache[sm] = f
inp = inp[~inp["smiles"].isin(bad_smiles)].reset_index(drop=True)
print(f"compounds featurised: {len(comp_cache)} ok, {len(bad_smiles)} unparseable; "
      f"rows now {len(inp)} ({inp['rxn_id'].nunique()} reactions)", flush=True)
prot_cache = {s: protein_words(s) for s in inp["seq"].unique()}
print(f"proteins featurised: {len(prot_cache)}", flush=True)

inp["Temp_K"] = inp["Temp_C"] + 273.15
inp["Temp_K_norm"] = (inp["Temp_K"] - TK_MIN) / (TK_MAX - TK_MIN)
inp["Inv_Temp_norm"] = (1.0 / inp["Temp_K"] - INV_MIN) / (INV_MAX - INV_MIN)

# ---------- model ----------
param = load_pkl(os.path.join(DLT, "data", "hyparams", "param_2.pkl"))  # latent=40, out=4: matches checkpoint
device = torch.device("cpu")
M = DLTKcat(len(fp_dict), len(word_dict), param["comp_dim"], param["prot_dim"],
            param["gat_dim"], param["num_head"], param["dropout"], param["alpha"],
            param["window"], param["layer_cnn"], param["latent_dim"], param["layer_out"])
M.to(device)
M.load_state_dict(torch.load(MODEL, map_location=device))
M.eval()
print("model loaded", flush=True)

def batch_pad(arr):
    N = max(a.shape[0] for a in arr)
    if arr[0].ndim == 1:
        out = np.zeros((len(arr), N)); mask = np.zeros((len(arr), N))
        for i, a in enumerate(arr):
            out[i, :a.shape[0]] = a; mask[i, :a.shape[0]] = 1
    else:
        out = np.zeros((len(arr), N, N)); mask = np.zeros((len(arr), N, N))
        for i, a in enumerate(arr):
            n = a.shape[0]; out[i, :n, :n] = a; mask[i, :n, :n] = 1
    return out, mask

BATCH = 64
preds = np.zeros(len(inp))
order = np.argsort([comp_cache[s][0].shape[0] for s in inp["smiles"]])  # size-sorted batches
with torch.no_grad():
    for bi_ in range(math.ceil(len(order) / BATCH)):
        idx = order[bi_ * BATCH:(bi_ + 1) * BATCH]
        rows = inp.iloc[idx]
        comps = [comp_cache[s] for s in rows["smiles"]]
        atoms_pad, atoms_mask = batch_pad([c[0] for c in comps])
        adj_pad, _ = batch_pad([c[1] for c in comps])
        fps = np.stack([c[2] for c in comps]).astype(float)
        amino_pad, amino_mask = batch_pad([prot_cache[s] for s in rows["seq"]])
        t = lambda x, ty: ty(x).to(device)
        pred = M(t(atoms_pad, torch.LongTensor), t(atoms_mask, torch.FloatTensor),
                 t(adj_pad, torch.LongTensor), t(amino_pad, torch.LongTensor),
                 t(amino_mask, torch.FloatTensor), t(fps, torch.FloatTensor),
                 t(rows[["Inv_Temp_norm"]].values, torch.FloatTensor),
                 t(rows[["Temp_K_norm"]].values, torch.FloatTensor))
        preds[idx] = pred.cpu().numpy().reshape(-1)
        if bi_ % 50 == 0:
            print(f"batch {bi_}/{math.ceil(len(order)/BATCH)}", flush=True)

inp["pred_log10kcat"] = preds
inp[["rxn_id", "enz", "sub", "Temp_C", "pred_log10kcat"]].to_csv(OUT, index=False)
print(f"wrote {OUT}: {len(inp)} predictions, {inp['rxn_id'].nunique()} reactions, "
      f"log10kcat range {preds.min():.2f}..{preds.max():.2f} "
      f"(median {np.median(preds):.2f})", flush=True)
