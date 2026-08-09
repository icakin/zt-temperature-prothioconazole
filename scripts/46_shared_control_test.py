#!/usr/bin/env python3
# Does the Fig. 5B slope pattern survive removal of shared-control covariance?
# Warming contrast (27C vs 15C control) and drug contrast (2 vs 0 at each T)
# share control samples at 15 C and 27 C but not at 21 C. Shared noise inflates
# the 15 C slope and deflates the 27 C slope. Test with disjoint control splits.
import numpy as np, pandas as pd
rng = np.random.default_rng(1)
V = pd.read_csv('tables/revision/expression/vsd_matrix.csv', index_col=0)
C = pd.read_csv('tables/revision/expression/coldata.csv')
C['temperature']=C['temperature'].astype(int); C['prothioconazole']=C['prothioconazole'].astype(int)
g = lambda T,d: C[(C.temperature==T)&(C.prothioconazole==d)]['sample'].tolist()
def m(cols): return V[cols].mean(axis=1)

# --- naive (shared-control) slopes, as in Fig 5B ---
warm = m(g(27,0)) - m(g(15,0))
naive = {}
for T in (15,21,27):
    drug = m(g(T,2)) - m(g(T,0))
    naive[T] = np.polyfit(warm, drug, 1)[0]
print("naive (shared controls):", {k: round(v,3) for k,v in naive.items()})

# --- disjoint-control splits ---
res = {15:[], 21:[], 27:[]}
NS = 200
for it in range(NS):
    c15 = g(15,0); c27 = g(27,0)
    rng.shuffle(c15); rng.shuffle(c27)
    a15,b15 = c15[:len(c15)//2], c15[len(c15)//2:]
    a27,b27 = c27[:len(c27)//2], c27[len(c27)//2:]
    warm_s = m(a27) - m(a15)                 # warming uses subset A
    for T in (15,21,27):
        ctrl = b15 if T==15 else (b27 if T==27 else g(21,0))   # drug uses subset B
        drug_s = m(g(T,2)) - m(ctrl)
        res[T].append(np.polyfit(warm_s, drug_s, 1)[0])
print("\ndisjoint-control splits (%d iterations):" % NS)
out=[]
for T in (15,21,27):
    a=np.array(res[T]); lo,md,hi=np.percentile(a,[2.5,50,97.5])
    print(f"  {T} C: median slope {md:.3f}   95%% range [{lo:.3f}, {hi:.3f}]   (naive {naive[T]:.3f})")
    out.append({"temperature":T,"slope_naive":naive[T],"slope_split_median":md,
                "slope_split_lo":lo,"slope_split_hi":hi})
d = np.array(res[15])-np.array(res[27])
print(f"\n  slope(15 C) - slope(27 C) under disjoint controls: median {np.median(d):.3f} "
      f"95% [{np.percentile(d,2.5):.3f}, {np.percentile(d,97.5):.3f}]")
print(f"  fraction of splits with slope(15) > slope(27): {100*np.mean(np.array(res[15])>np.array(res[27])):.1f}%")
print(f"  fraction with slope(21) > slope(27):          {100*np.mean(np.array(res[21])>np.array(res[27])):.1f}%")
pd.DataFrame(out).to_csv('tables/gem/shared_control_test.csv', index=False)
