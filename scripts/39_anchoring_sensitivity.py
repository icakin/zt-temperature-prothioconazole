#!/usr/bin/env python3
# 39_anchoring_sensitivity.py — emergent TPC with enzyme-optimum anchors shifted.
# Requires the etcGEMs package installed (pip install -e etcGEMs) and the
# zt_ipo323 strain folder. Run from the project root.
import pandas as pd, subprocess, shutil, os
TH = "etcGEMs/strains/zt_ipo323/thermal/BestParamsTopt.csv"
OUTD = "etcGEMs/strains/zt_ipo323/outputs/tpc/nominal_tpc.csv"
orig = pd.read_csv(TH)
shutil.copy(TH, TH + ".orig")
try:
    for sh, tag in [(-4, "-4"), (-2, "-2"), (0, "0"), (+2, "+2")]:
        b = orig.copy(); b["Topt"] = b["Topt"] + sh
        b.to_csv(TH, index=False)
        subprocess.run(["etcgem", "tpc", "--strain", "zt_ipo323"], check=True)
        shutil.copy(OUTD, f"tables/gem/anchor_tpc_shift{tag}.csv")
        print("anchor shift", tag, "done")
finally:
    shutil.copy(TH + ".orig", TH)
    subprocess.run(["etcgem", "tpc", "--strain", "zt_ipo323"], check=True)
print("originals restored")
