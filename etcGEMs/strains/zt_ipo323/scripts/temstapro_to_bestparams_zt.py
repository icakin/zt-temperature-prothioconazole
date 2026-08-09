#!/usr/bin/env python3
"""temstapro_to_bestparams_zt.py — turn TemStaPro --mean-output into
thermal/BestParamsTopt.csv for the zt_ipo323 etc-GEM strain.

Differences from the cauris version: Z. tritici has no homology meltome to
calibrate against, so the sequence pseudo-Tm is used directly (flagged in
TmTag). Topt is set as CENTER + SLOPE*(Tm - meanTm), same construction as
build_bestparams_topt.py.

Usage:
  python temstapro_to_bestparams_zt.py \
      --tsv thermal/enzyme_tm_temstapro_raw.tsv \
      --faa thermal/enzyme_sequences.faa \
      --out thermal/BestParamsTopt.csv \
      [--center 24 --slope 0.5 --dcp -4000]
Then in strain.yaml: uncomment enzyme_params and set thermal_model: unfolding.
"""
import argparse, numpy as np, pandas as pd

THRESH = np.array([40, 45, 50, 55, 60, 65], float)

def pseudo_tm(row):
    p = np.array([row[f"t{t:.0f}_raw"] for t in THRESH], float)
    if p[0] < 0.5:   return float(35.0 + 5.0*(p[0]/0.5))
    if p[-1] >= 0.5: return float(65.0 + 5.0*min((p[-1]-0.5)/0.5, 1.0))
    i = np.where(p >= 0.5)[0][-1]
    x0, x1, y0, y1 = THRESH[i], THRESH[i+1], p[i], p[i+1]
    return float(x0 + (x1-x0)*(y0-0.5)/(y0-y1)) if y0 != y1 else float(x0)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--center", type=float, default=24.0,
                    help="proteome-mean Topt in C (Z. tritici growth optimum ~24)")
    ap.add_argument("--slope", type=float, default=0.5,
                    help="Topt heterogeneity per degree of Tm (catalytic optima vary less than Tm)")
    ap.add_argument("--dcp", type=float, default=-4000.0,
                    help="MMRT dCpt prior, J/mol/K (literature)")
    a = ap.parse_args()

    df = pd.read_csv(a.tsv, sep="\t")
    idcol = df.columns[0]
    df["Tm_C"] = df.apply(pseudo_tm, axis=1)

    lengths = {}
    gid = None; n = 0
    for line in open(a.faa):
        if line.startswith(">"):
            if gid: lengths[gid] = n
            gid = line[1:].split()[0]; n = 0
        else:
            n += len(line.strip())
    if gid: lengths[gid] = n

    tm_mean = df["Tm_C"].mean()
    out = pd.DataFrame({
        "Topt": 273.15 + a.center + a.slope*(df["Tm_C"] - tm_mean),
        "Topt_std": 10.0,
        "Length": df[idcol].map(lengths).fillna(400).astype(int),
        "Tm": 273.15 + df["Tm_C"],
        "Tm_std": 5.0,
        "T90": np.nan,
        "dCpt": a.dcp,
        "dCpt_std": 1000.0,
        "topt_source": "temstapro_pseudoTm",
        "TmTag": "TemStaPro",
    })
    out.index = df[idcol].astype(str).str.replace(r"\s.*", "", regex=True)
    out.to_csv(a.out)
    print(f"wrote {a.out}: {len(out)} enzymes, "
          f"Tm {df['Tm_C'].min():.1f}-{df['Tm_C'].max():.1f} C (mean {tm_mean:.1f}), "
          f"Topt {out['Topt'].min()-273.15:.1f}-{out['Topt'].max()-273.15:.1f} C")

if __name__ == "__main__":
    main()
