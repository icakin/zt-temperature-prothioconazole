#!/usr/bin/env python3
# =============================================================================
# 31_ribosome_forecast.py — Prereg v4, Approach B (FROZEN).
# Control-trained ribosome forecast:
#   fit log(growth) = a + b*s on the 3 CONTROL conditions (s = ribosome module
#   score), predict the 3 drug conditions blind. Gates B1 (predictions inside
#   measured 95% CrI) and B2 (predicted attenuation d27 > d15).
# Secondary: leave-one-temperature-pair-out (fit on 4, predict held-out pair).
# Uncertainty: RNA-replicate bootstrap x growth Bayesian bootstrap (500 draws).
# =============================================================================
import os, re, numpy as np, pandas as pd
rng=np.random.default_rng(42)
OUT="tables/revision/v4"; os.makedirs(OUT,exist_ok=True)
LOG=open(f"{OUT}/84_log.txt","w")
def log(*a):
    s=" ".join(str(x) for x in a); print(s); LOG.write(s+"\n"); LOG.flush()

# ---- ribosome module score (identical definition to script 81) ---------------
cnt=pd.read_csv("tables/rnaseq_stage/normalized_counts.csv",index_col=0)
logn=np.log2(cnt+1)
de21=pd.read_csv("tables/rnaseq_stage/DE_proth_at_21C_2vs0.csv")
desc=de21["description"].fillna("").str.lower(); nm=de21["gene_name"].fillna("").str.lower()
kegg=de21["KEGG"].fillna("").astype(str); pfam=de21["PFAM"].fillna("").str.lower()
# ribosome module = KEGG ribosome pathway (map03010), identical to Fig 5C / script 81
# (prereg §7 amendment 1: module unified to map03010; verdict reported under both
# the map03010 set and the broader annotation-OR set to preclude tuning)
ribo=[g for g in de21[kegg.str.contains("map03010")]["gene_id"] if g in cnt.index]
log(f"ribosome module: {len(ribo)} genes")
def cols(T,D): return [c for c in cnt.columns if re.match(fr"12129_{T}_{D}_R",c)]
CONDS=[("15","0"),("21","0"),("27","0"),("15","2"),("21","2"),("27","2")]
score_samples={(T,D):logn.loc[ribo,cols(T,D)].mean(axis=0) for T,D in CONDS}
score={k:float(v.mean()) for k,v in score_samples.items()}

# ---- measured growth: replicate vials at the six cells ------------------------
# (clarification logged in prereg §7: measured-growth distribution = Bayesian
# bootstrap over replicate vials at the exact T x dose cells, the same source
# used by the condition-constrained flux analysis)
d=pd.read_csv("tables/physiology/derived_N0_R_results_with_carbon.csv")
d["conc"]=d["Dose"].replace("Control","0").astype(float); d["Tc"]=d["T"].astype(int)
d=d[(d.growth_C_per_C_h>0)&(d.respiration_C_per_C_h>0)]
gobs={(T,D):d[(d.Tc==int(T))&(d.conc==float(D))].growth_C_per_C_h.values for T,D in CONDS}
for k,v in gobs.items(): log(f"growth obs {k}: n={len(v)} mean={v.mean():.4f}")

def bb_mean(v,n=4000):
    """Bayesian bootstrap distribution of the mean."""
    w=rng.dirichlet(np.ones(len(v)),size=n)
    return w@v
gdist={k:bb_mean(v) for k,v in gobs.items()}
gmed={k:float(np.median(vd)) for k,vd in gdist.items()}
gci={k:(float(np.percentile(vd,2.5)),float(np.percentile(vd,97.5))) for k,vd in gdist.items()}

# ---- PRIMARY: fit on 3 controls, predict 3 drug conditions --------------------
CTRL=[("15","0"),("21","0"),("27","0")]; DRUG=[("15","2"),("21","2"),("27","2")]
def fit_predict(sc,gm,train,test):
    x=np.array([sc[k] for k in train]); y=np.log([gm[k] for k in train])
    b,a=np.polyfit(x,y,1)
    return {k:float(np.exp(a+b*sc[k])) for k in test},float(b),float(a)
pred,b,a=fit_predict(score,gmed,CTRL,DRUG)
log(f"\ncontrol fit: log(g) = {a:.3f} + {b:.3f}*s  (n=3)")
rows=[]
for k in DRUG:
    lo,hi=gci[k]; inside=lo<=pred[k]<=hi
    rows.append(dict(T=k[0],D=k[1],pred=pred[k],obs=gmed[k],lo95=lo,hi95=hi,
                     inside=inside,resid_log=float(np.log(pred[k]/gmed[k]))))
    log(f"predict {k}: pred={pred[k]:.4f}  obs={gmed[k]:.4f} [{lo:.4f},{hi:.4f}]  inside={inside}  log-resid={np.log(pred[k]/gmed[k]):+.3f}")
B1=all(r["inside"] for r in rows)
d15=np.log(pred[("15","2")]/gmed[("15","0")]); d27=np.log(pred[("27","2")]/gmed[("27","0")])
B2=bool(d27>d15)
log(f"predicted delta15={d15:+.3f}  delta27={d27:+.3f}  -> B2 (d27>d15): {B2}")
log(f"GATE B1 (all 3 inside 95% CrI): {B1}")
log(f"GATE B2 (attenuation reproduced): {B2}")

# ---- gate robustness: joint RNA-bootstrap x growth-bootstrap ------------------
B=500; b1h=0; b2h=0
for i in range(B):
    sc_b={}
    for k in CONDS:
        ss=score_samples[k].values; sc_b[k]=float(rng.choice(ss,len(ss),replace=True).mean())
    gm_b={k:float(rng.choice(gobs[k],len(gobs[k]),replace=True).mean()) for k in CONDS}
    try: pr,_,_=fit_predict(sc_b,gm_b,CTRL,DRUG)
    except Exception: continue
    ok1=all(gci[k][0]<=pr[k]<=gci[k][1] for k in DRUG)
    dd15=np.log(pr[("15","2")]/gm_b[("15","0")]); dd27=np.log(pr[("27","2")]/gm_b[("27","0")])
    b1h+=ok1; b2h+=(dd27>dd15)
log(f"robustness (500 joint bootstraps): P(B1)={b1h/B:.3f}  P(B2)={b2h/B:.3f}")

# ---- SECONDARY: leave-one-temperature-pair-out --------------------------------
log("\nleave-one-temperature-pair-out (fit on 4 conditions, predict the held-out pair):")
rows2=[]
for Thold in ["15","21","27"]:
    train=[k for k in CONDS if k[0]!=Thold]; test=[k for k in CONDS if k[0]==Thold]
    pr,bb,aa=fit_predict(score,gmed,train,test)
    for k in test:
        lo,hi=gci[k]
        rows2.append(dict(hold=Thold,T=k[0],D=k[1],pred=pr[k],obs=gmed[k],lo95=lo,hi95=hi,
                          inside=lo<=pr[k]<=hi,resid_log=float(np.log(pr[k]/gmed[k])),slope=bb))
        log(f"  hold {Thold}: {k} pred={pr[k]:.4f} obs={gmed[k]:.4f} [{lo:.4f},{hi:.4f}] inside={lo<=pr[k]<=hi} log-resid={np.log(pr[k]/gmed[k]):+.3f}")
# ribosome-only LOTPO total squared log-residual (baseline for gate A4)
ssr=float(np.sum([r["resid_log"]**2 for r in rows2]))
log(f"ribosome-only LOTPO total squared log-residual (A4 baseline): {ssr:.4f}")

pd.DataFrame(rows).to_csv(f"{OUT}/B_primary_predictions.csv",index=False)
pd.DataFrame(rows2).to_csv(f"{OUT}/B_lotpo_predictions.csv",index=False)
pd.DataFrame(dict(cond=[f"{t}_{dd}" for t,dd in CONDS],score=[score[k] for k in CONDS],
                  growth_med=[gmed[k] for k in CONDS],
                  growth_lo=[gci[k][0] for k in CONDS],growth_hi=[gci[k][1] for k in CONDS]))\
  .to_csv(f"{OUT}/B_condition_summary.csv",index=False)
log(f"\nVERDICT B: {'PASS' if (B1 and B2) else 'FAIL'} (B1={B1}, B2={B2})")
LOG.close()
