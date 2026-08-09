"""Augment the tuned control / identifiability tables with local enzyme identities.

Adds a gene symbol, a short enzyme/protein name and (where available) an EC number to
each rxn_id/enzyme_id row, sourced ENTIRELY from local files — no web fetch:

* enzyme/protein NAME  <- the ecModel reaction .name (GECKO isozyme suffixes stripped);
* gene SYMBOL          <- the measured proteomics table (b-number -> genename), joined
                          through the model's reaction->gene (b-number) map;
* UniProt accession     <- the enzyme_id already in the table;
* EC number             <- the ecModel reaction annotation if present (it is not stored
                          in this GECKO build, so left blank rather than invented).

Writes *_annotated.csv next to the source tables. Idempotent; re-run any time.
"""
from __future__ import annotations

import os
import re
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
MODEL = os.path.join(ROOT, "strains", "eciML1515", "model", "eciML1515_batch.xml")
PROT = os.path.join(ROOT, "strains", "eciML1515", "proteomics", "tem_proteomic.csv")
CTRL = os.path.join(ROOT, "strains", "eciML1515", "outputs", "control_tuned")


def _clean_name(name: str) -> str:
    s = re.sub(r"\s*\(No\d+\)\s*$", "", str(name))          # drop GECKO isozyme tag
    s = s.replace(" (reversible)", "")
    return s.strip()


def _rxn_base(rid: str) -> str:
    return re.sub(r"(No\d+)?(_REV)?$", "", str(rid))


def build_maps():
    sys.path.insert(0, os.path.join(ROOT, "src"))
    from etcgem.providers import _load_any
    m = _load_any(MODEL)
    prot = pd.read_csv(PROT)
    bgene = {a: str(g) for a, g in zip(prot["Accession"], prot["genename"]) if isinstance(g, str)}

    def annot(rid):
        r = None
        for cand in (rid, _rxn_base(rid)):
            if cand in m.reactions:
                r = m.reactions.get_by_id(cand)
                break
        if r is None:
            return "", "", ""
        bs = [g.id for g in r.genes]
        gene = ";".join(sorted({bgene[b] for b in bs if b in bgene}))
        ec = r.annotation.get("ec-code", "")
        if isinstance(ec, list):
            ec = ";".join(ec)
        return gene, _clean_name(r.name), (ec or "")
    return annot


def augment(path_in, path_out, annot, front_cols):
    df = pd.read_csv(path_in)
    ann = df["rxn_id"].map(lambda r: annot(r))
    df.insert(df.columns.get_loc("rxn_id") + 1, "gene", [a[0] for a in ann])
    df.insert(df.columns.get_loc("gene") + 1, "enzyme_name", [a[1] for a in ann])
    df["EC"] = [a[2] for a in ann]
    # bring the identity columns to the front for readability
    cols = [c for c in front_cols if c in df.columns] + \
           [c for c in df.columns if c not in front_cols]
    df[cols].to_csv(path_out, index=False)
    return len(df), int((df["gene"] != "").sum())


ENV_CC = ["CC[Topt_C,Topt_i]", "CC[Topt_C,dCp_i]", "CC[CT_max_C,Topt_i]",
          "CC[CT_max_C,dCp_i]", "CC[niche_width_C,Topt_i]", "CC[niche_width_C,dCp_i]"]


def _short_name(name: str) -> str:
    # drop a trailing substrate parenthetical for a cleaner display label
    return re.sub(r"\s*\([^)]*\)\s*(\(r\))?\s*$", "", str(name)).strip() or str(name)


def rerank_and_plot():
    """Add the ENVELOPE CONTROL COEFFICIENT (max |CC| over the envelope descriptors x
    {Topt_i, dCp_i}) to the annotated control table, re-rank by it, and regenerate the
    'which enzymes matter' figure + a deduped top-enzyme table (both from the SAME data
    so figure and table are consistent)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    p = os.path.join(CTRL, "thermal_control_annotated.csv")
    d = pd.read_csv(p)
    cc = [c for c in ENV_CC if c in d.columns]
    d["control_coeff"] = d[cc].abs().max(axis=1)
    d = d.sort_values("control_coeff", ascending=False).reset_index(drop=True)
    d["rank"] = range(1, len(d) + 1)
    cols = ["rank", "gene", "enzyme_name", "enzyme_id", "rxn_id", "EC",
            "control_coeff", "thermal_screen"] + cc
    d[[c for c in cols if c in d.columns]].to_csv(p, index=False)

    # deduped top enzymes (one row per enzyme_id = its strongest reaction), for display
    top = (d.sort_values("control_coeff", ascending=False)
             .drop_duplicates("enzyme_id").head(15).reset_index(drop=True))
    top.insert(0, "display_rank", range(1, len(top) + 1))
    top["enzyme_short"] = top["enzyme_name"].map(_short_name)

    def _label(r):
        g = str(r["gene"]).strip()
        if g and g.lower() != "nan":
            return g
        en = str(r["enzyme_short"]).strip()
        return en if en and en.lower() != "nan" else str(r["rxn_id"])
    top["label"] = top.apply(_label, axis=1)
    top_path = os.path.join(CTRL, "control_top_enzymes.csv")
    top[["display_rank", "gene", "enzyme_short", "enzyme_id", "control_coeff",
         "CC[niche_width_C,Topt_i]", "CC[CT_max_C,Topt_i]"]].rename(
        columns={"enzyme_short": "enzyme_name"}).to_csv(top_path, index=False)

    # figure: single clean control-coefficient series, gene-labelled
    fig, ax = plt.subplots(figsize=(9, 4.6))
    y = top["control_coeff"].values
    x = range(len(top))
    ax.bar(x, y, color="tab:red", alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels(top["label"].astype(str).str.slice(0, 14), rotation=55, ha="right", fontsize=8)
    ax.set_ylabel("envelope control coefficient\n(max |∂descriptor/∂param|)")
    ax.set_title("Which enzymes control the thermal envelope (tuned model, rich BHI)\n"
                 "carried by the niche-width response; per-enzyme control of $T_{opt}$ is ~0", fontsize=10)
    fig.tight_layout()
    out = os.path.join(CTRL, "control_coefficient_bar.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"[annotate] re-ranked by control_coeff; top enzyme {top.iloc[0]['label']} "
          f"(cc={top.iloc[0]['control_coeff']:.3f}); wrote {os.path.basename(top_path)} + figure")


def main():
    annot = build_maps()
    n1, g1 = augment(os.path.join(CTRL, "thermal_control.csv"),
                     os.path.join(CTRL, "thermal_control_annotated.csv"), annot,
                     ["rank", "gene", "enzyme_name", "enzyme_id", "rxn_id", "EC", "thermal_screen"])
    n2, g2 = augment(os.path.join(CTRL, "identifiability.csv"),
                     os.path.join(CTRL, "identifiability_annotated.csv"), annot,
                     ["gene", "enzyme_name", "enzyme_id", "rxn_id", "EC", "parameter",
                      "ident", "identifiable_from_growth", "refined"])
    print(f"[annotate] thermal_control: {n1} rows, {g1} with a gene symbol")
    print(f"[annotate] identifiability: {n2} rows, {g2} with a gene symbol")
    rerank_and_plot()


if __name__ == "__main__":
    main()
