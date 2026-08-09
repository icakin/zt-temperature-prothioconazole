#!/usr/bin/env python3
# =============================================================================
# 49_climate_relevance.py — how much of the European wheat fungicide-spray
# window now falls in the thermal bands where we measured antagonism?
#
# The physiology in this study is bounded: we assayed 15-28 C, found synergy
# below the thermal optimum and antagonism above it. This script asks the only
# question that turns that into an agronomic statement: has the number of
# spray-window hours spent in each of those measured bands changed?
#
# Data   : ERA5 reanalysis hourly 2 m air temperature, served by the Open-Meteo
#          historical archive API (free, no API key, CC-BY-4.0).
#          https://open-meteo.com/en/docs/historical-weather-api
# Sites  : ten representative NW-European winter-wheat regions where Septoria
#          tritici blotch is economically important.
# Window : 1 April - 15 July, spanning the T1-T3 fungicide timings (GS31-GS59)
#          on winter wheat in maritime Europe.
# Periods: 1961-1990 (baseline normal) vs 1996-2025 (recent normal).
#
# IMPORTANT CAVEATS, carried into the manuscript text:
#   * 2 m air temperature is not leaf or canopy temperature; canopy conditions
#     inside a closed wheat crop are typically cooler and more humid by day.
#   * Hours are not weighted by infection risk or by when sprays are actually
#     applied; this is exposure accounting, not a disease model.
#   * Bands are the empirically assayed range only. Hours below 15 C or above
#     28 C are reported separately and are NOT interpreted.
#
# RUN FROM the master "X0123 copy" folder:
#     python3 scripts/49_climate_relevance.py
# Requires network access. First run downloads ~600 small JSON responses
# (2-5 min) and caches them under data/climate/cache/; re-runs are instant.
# Outputs -> tables/revision/climate/
# =============================================================================
import os, sys, json, time, math
import urllib.request, urllib.error, urllib.parse
import numpy as np
import pandas as pd

CACHE = "data/climate/cache"
OUT = "tables/revision/climate"
os.makedirs(CACHE, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

API = "https://archive-api.open-meteo.com/v1/archive"

# --- sites: NW-European winter-wheat regions with high STB pressure ----------
SITES = [
    ("Rothamsted, UK",        51.81,  -0.36),
    ("Cambridge, UK",         52.20,   0.13),
    ("Oak Park, IE",          52.86,  -6.92),
    ("Rennes, FR",            48.11,  -1.68),
    ("Orleans (Beauce), FR",  47.90,   1.90),
    ("Ghent, BE",             51.05,   3.72),
    ("Wageningen, NL",        51.97,   5.67),
    ("Halle, DE",             51.50,  11.97),
    ("Roskilde, DK",          55.65,  12.09),
    ("Uppsala, SE",           59.86,  17.64),
]

# --- spray window and periods ------------------------------------------------
WIN_START, WIN_END = "04-01", "07-15"
EARLY = (1961, 1990)
RECENT = (1996, 2025)
YEARS = list(range(EARLY[0], EARLY[1] + 1)) + list(range(RECENT[0], RECENT[1] + 1))

# --- thermal bands, anchored on what was actually measured -------------------
# 24.0 C is the measured thermal optimum and the Bliss zero-deviation reference.
BANDS = [
    ("below_assayed", -math.inf, 15.0),   # outside the assayed range
    ("synergy_15_20",       15.0, 20.0),  # measured Bliss deviation < 0
    ("neutral_20_24",       20.0, 24.0),  # transitional, straddles T_opt
    ("antagonism_24_28",    24.0, 28.0),  # measured Bliss deviation > 0
    ("above_assayed",       28.0, math.inf),
]
BAND_NAMES = [b[0] for b in BANDS]


def fetch_year(lat, lon, year, tries=4):
    """Hourly 2 m temperature for one site-year over the spray window."""
    key = f"{lat:.2f}_{lon:.2f}_{year}.json"
    path = os.path.join(CACHE, key)
    if os.path.exists(path):
        with open(path) as fh:
            return json.load(fh)
    q = urllib.parse.urlencode({
        "latitude": lat, "longitude": lon,
        "start_date": f"{year}-{WIN_START}", "end_date": f"{year}-{WIN_END}",
        "hourly": "temperature_2m", "timezone": "UTC",
    })
    url = f"{API}?{q}"
    for attempt in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=90) as r:
                js = json.load(r)
            if "hourly" not in js:
                raise RuntimeError(js.get("reason", "no hourly block"))
            with open(path, "w") as fh:
                json.dump(js, fh)
            time.sleep(0.2)              # be polite to the free API
            return js
        except Exception as e:                       # rate limit / transient
            if attempt == tries - 1:
                raise
            time.sleep(2.0 * (attempt + 1))
    raise RuntimeError("unreachable")


def band_counts(temps):
    t = np.asarray(temps, dtype=float)
    t = t[np.isfinite(t)]
    out = {"n_hours": int(t.size), "mean_C": float(t.mean()) if t.size else np.nan}
    for name, lo, hi in BANDS:
        out[name] = int(((t >= lo) & (t < hi)).sum())
    return out


# --- download and tally ------------------------------------------------------
rows = []
total = len(SITES) * len(YEARS)
done = 0
for site, lat, lon in SITES:
    for yr in YEARS:
        js = fetch_year(lat, lon, yr)
        c = band_counts(js["hourly"]["temperature_2m"])
        c.update(site=site, lat=lat, lon=lon, year=yr,
                 period="early" if yr <= EARLY[1] else "recent")
        rows.append(c)
        done += 1
        if done % 50 == 0:
            print(f"  {done}/{total} site-years", flush=True)
S = pd.DataFrame(rows)
S.to_csv(os.path.join(OUT, "climate_site_year_bands.csv"), index=False)
print(f"site-years: {len(S)}  ({S.year.min()}-{S.year.max()}, {S.site.nunique()} sites)")

# --- per-site early vs recent ------------------------------------------------
per_site = (S.groupby(["site", "period"])[BAND_NAMES + ["mean_C", "n_hours"]]
              .mean().reset_index())
piv = per_site.pivot(index="site", columns="period")
recs = []
for site in S.site.unique():
    r = {"site": site}
    for b in BAND_NAMES + ["mean_C"]:
        e, n = piv[(b, "early")][site], piv[(b, "recent")][site]
        r[f"{b}_early"] = round(e, 1)
        r[f"{b}_recent"] = round(n, 1)
        r[f"{b}_change"] = round(n - e, 1)
    recs.append(r)
P = pd.DataFrame(recs)
P.to_csv(os.path.join(OUT, "climate_per_site_change.csv"), index=False)

# --- pooled summary, with a paired test across sites -------------------------
from scipy import stats
pool = []
for b in BAND_NAMES + ["mean_C"]:
    e = P[f"{b}_early"].values
    n = P[f"{b}_recent"].values
    t = stats.ttest_rel(n, e)
    pool.append({
        "band": b,
        "early_mean_hours": round(e.mean(), 1),
        "recent_mean_hours": round(n.mean(), 1),
        "change_hours": round(n.mean() - e.mean(), 1),
        "pct_change": round(100 * (n.mean() - e.mean()) / e.mean(), 1) if e.mean() else np.nan,
        "sites_increasing": int((n > e).sum()),
        "n_sites": len(e),
        "paired_t": round(float(t.statistic), 3),
        "paired_p": float(t.pvalue),
    })
POOL = pd.DataFrame(pool)
POOL.to_csv(os.path.join(OUT, "climate_pooled_summary.csv"), index=False)

# --- linear trend in antagonism-band hours, per site and pooled --------------
tr = []
for site, g in S.groupby("site"):
    sl = stats.linregress(g["year"], g["antagonism_24_28"] + g["above_assayed"])
    tr.append({"site": site, "slope_hours_per_decade": round(sl.slope * 10, 2),
               "p": float(sl.pvalue), "r2": round(sl.rvalue ** 2, 3)})
ann = S.groupby("year")[["antagonism_24_28", "above_assayed",
                         "synergy_15_20", "mean_C"]].mean().reset_index()
ann["supra_optimal"] = ann["antagonism_24_28"] + ann["above_assayed"]
slp = stats.linregress(ann["year"], ann["supra_optimal"])
tr.append({"site": "ALL SITES (mean)", "slope_hours_per_decade": round(slp.slope * 10, 2),
           "p": float(slp.pvalue), "r2": round(slp.rvalue ** 2, 3)})
TR = pd.DataFrame(tr)
TR.to_csv(os.path.join(OUT, "climate_trend_supraoptimal.csv"), index=False)
ann.to_csv(os.path.join(OUT, "climate_annual_pooled.csv"), index=False)

# --- report ------------------------------------------------------------------
lines = []
w = lines.append
w("# Climate relevance of the measured thermal bands\n")
w(f"ERA5 hourly 2 m air temperature, {len(SITES)} NW-European winter-wheat sites,")
w(f"spray window {WIN_START} to {WIN_END}, {EARLY[0]}-{EARLY[1]} vs {RECENT[0]}-{RECENT[1]}.\n")
w("## Mean spray-window hours per season, per site\n")
w(POOL.to_string(index=False))
w("")
sup_e = POOL.loc[POOL.band == "antagonism_24_28", "early_mean_hours"].iat[0] + \
    POOL.loc[POOL.band == "above_assayed", "early_mean_hours"].iat[0]
sup_r = POOL.loc[POOL.band == "antagonism_24_28", "recent_mean_hours"].iat[0] + \
    POOL.loc[POOL.band == "above_assayed", "recent_mean_hours"].iat[0]
w(f"\nSupra-optimal (>24 C) spray-window hours: {sup_e:.0f} -> {sup_r:.0f} per season "
  f"({100*(sup_r-sup_e)/sup_e:+.0f}%).")
w(f"Pooled trend: {slp.slope*10:+.1f} h per decade (p = {slp.pvalue:.2g}, "
  f"R2 = {slp.rvalue**2:.2f}).\n")
w("## Per-site trend in supra-optimal hours\n")
w(TR.to_string(index=False))
w("\n## Caveats\n")
w("2 m air temperature is not canopy temperature. Hours are unweighted by")
w("infection risk or actual spray timing. Bands outside 15-28 C are reported")
w("but not interpreted, as they were not assayed.")
rep = "\n".join(lines)
with open(os.path.join(OUT, "climate_relevance_report.md"), "w") as fh:
    fh.write(rep + "\n")
print("\n" + rep)
print(f"\nwrote -> {OUT}")
