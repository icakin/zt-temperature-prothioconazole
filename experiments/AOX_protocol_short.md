# AOX-inhibitor × prothioconazole — short protocol
*Z. tritici* IPO323 · PreSens SensorDish · 5 mL vials · 0.4% DMSO, 0.2% ethanol

**Question:** does nPG (confirmed by SHAM) lower the prothioconazole EC₅₀ at 27 °C
but not at 15 °C? Signature = inhibitor × temperature interaction on EC₅₀.

**Plates:** 3 replicates × (one 15 °C + one 27 °C plate) = 6 growth plates
+ 1 azide-validation plate. Both temperatures come from the SAME master, same day.

---

## 1. Stocks

Prothioconazole ladder = serial 1:1 dilutions of your **1 mg/mL** stock in DMSO,
added at **20 µL/vial** (→ 0.4% DMSO, matches your assay):

| Stock (20 µL) | = mg/mL | → in vial |
|---|---|---|
| P4 | 1 (your stock, neat) | 4 mg/L |
| P2 | 0.5 | 2 mg/L |
| P1 | 0.25 | 1 mg/L |
| P0.5 | 0.125 | 0.5 mg/L |
| P0.25 | 0.0625 | 0.25 mg/L |
| P0 | neat DMSO | 0 mg/L |

| Inhibitor | Make | Solvent |
|---|---|---|
| nPG | 50 mM (21.2 mg → 2 mL) | ethanol |
| SHAM | 250 mM (76.6 mg → 2 mL) | ethanol |
| Azide | 200 mM (26 mg → 2 mL) | water |

Drug: **20 µL/vial** → 0.4% DMSO (no solubility worry — your 1 mg/mL stock is already
dissolved). Inhibitor: **10 µL-equivalent/vial** → 0.2% ethanol. Confirm SHAM (250 mM)
dissolves in ethanol; foil-wrap nPG/SHAM.

## 2. Precultures

3 independent precultures (R1–R3), mid-exponential. On the day, make **15 mL at
OD₆₀₀ 0.050** from each (µL preculture = 0.050 × 15000 ÷ your OD; top up with medium).
Use 2 mL for the growth master, 10 mL for the azide flask.

## 3. Per replicate — cups (do R1, R2, R3 separately)

| Cup | Recipe | Fills |
|---|---|---|
| **Master** (200 mL) | 197.995 mL medium + 2.005 mL OD-0.05 cells | → the arm cups |
| **Vehicle arm** (63 mL) | 62.87 mL master + 126 µL ethanol | 6 vials × 2 plates |
| **nPG arm** (63 mL) | 62.87 mL master + 126 µL 50 mM nPG | 6 × 2 |
| **SHAM arm** (63 mL) | 62.87 mL master + 126 µL 250 mM SHAM | 6 × 2 |
| **Vehicle blank** (21 mL) | 20.96 mL sterile medium + 42 µL ethanol | 2 × 2 |
| **nPG blank** (21 mL) | 20.96 mL sterile medium + 42 µL 50 mM nPG | 2 × 2 |
| **SHAM blank** (21 mL) | 20.96 mL sterile medium + 42 µL 250 mM SHAM | 2 × 2 |

## 4. Growth plate — what goes in each vial (same map for 15 °C and 27 °C)

Every vial = **4.98 mL cup mix + 20 µL drug stock.**

| Vial | Cup mix (4.98 mL) | Drug (20 µL) | → in vial |
|---|---|---|---|
| A1 | vehicle arm | P0 | 0 mg/L, no inhib |
| A2 | vehicle arm | P0.25 | 0.25 |
| A3 | vehicle arm | P0.5 | 0.5 |
| A4 | vehicle arm | P1 | 1 |
| A5 | vehicle arm | P2 | 2 |
| A6 | vehicle arm | P4 | 4 |
| B1 | nPG arm | P0 | 0, +nPG 0.1 mM |
| B2 | nPG arm | P0.25 | 0.25, +nPG |
| B3 | nPG arm | P0.5 | 0.5, +nPG |
| B4 | nPG arm | P1 | 1, +nPG |
| B5 | nPG arm | P2 | 2, +nPG |
| B6 | nPG arm | P4 | 4, +nPG |
| C1 | SHAM arm | P0 | 0, +SHAM 0.5 mM |
| C2 | SHAM arm | P0.25 | 0.25, +SHAM |
| C3 | SHAM arm | P0.5 | 0.5, +SHAM |
| C4 | SHAM arm | P1 | 1, +SHAM |
| C5 | SHAM arm | P2 | 2, +SHAM |
| C6 | SHAM arm | P4 | 4, +SHAM |
| D1 | vehicle blank | P0 | blank |
| D2 | vehicle blank | P0 | blank |
| D3 | nPG blank | P0 | blank +nPG |
| D4 | nPG blank | P0 | blank +nPG |
| D5 | SHAM blank | P0 | blank +SHAM |
| D6 | SHAM blank | P0 | blank +SHAM |

Rows A–C = cells; row D = cell-free (no cells). Fill the 15 °C and 27 °C plate from
the same cups, seal, start both readers together, run ~4 days (dark). All 3
replicates set up the same day.

## 5. Azide validation — separate, dense cells, acute (fills one 24-vial plate)

**Flask per replicate:** 10.000 mL OD-0.05 cells + ~39.8 mL medium + 200 µL P2
(0.5 mg/mL) = 50 mL → OD 0.01, 2 mg/L drug, 0.4% DMSO. Grow at 27 °C to
**OD₆₀₀ 0.15–0.30** (1–2 days later). Standardise R1–R3 to the same OD with medium
containing 2 mg/L drug + 0.4% DMSO. (Dense cells already carry the drug, so the
acute vials add no further drug.)

Every vial = **4.965 mL matrix + 10 µL inhibitor**, then inject **25 µL** after a
stable ~10 min baseline; record 15–20 min after. Final: azide 1 mM, nPG 0.1 mM,
SHAM 0.5 mM, ethanol 0.2%.

**8 vials per replicate** (R1 = cols 1–2, R2 = cols 3–4, R3 = cols 5–6):

| Row | Matrix (4.965 mL) | Inhibitor (10 µL) | odd col inject | even col inject |
|---|---|---|---|---|
| A | dense cells | ethanol | 25 µL water | 25 µL azide |
| B | dense cells | 50 mM nPG | 25 µL water | 25 µL azide |
| C | dense cells | 250 mM SHAM | 25 µL water | 25 µL azide |
| D | cell-free medium | nPG (odd) / SHAM (even) | 25 µL azide | 25 µL azide |

So per replicate: A = ±azide vehicle, B = ±azide nPG, C = ±azide SHAM (the water
column is the injection/inhibitor-only control), D = nPG blank + SHAM blank (both
azide). 8 × 3 = 24 vials, one plate.

**Read:** azide-resistant respiration = the residual after azide in the vehicle
vial (A even); nPG (B even) and SHAM (C even) should each cut it below that; the
water columns and the two D blanks separate injection artifact and chemical O₂ use.

## 6. Analysis

- Correct each growth trace with its matched blank (D1-2→A, D3-4→B, D5-6→C).
- Hill fit per temperature × arm → EC₅₀ (shared slope; the 27 °C control EC₅₀ may sit
  above 4 mg/L — report it as a wide/one-sided CI, not a point value).
- Model: log EC₅₀ ~ temperature × inhibitor + (1|replicate). Look for lower EC₅₀ with
  nPG **and** SHAM at 27 °C, little change at 15 °C, significant interaction.
- Azide: residual = post/pre-injection rate; nPG and SHAM should both cut the
  azide-resistant residual.
- **Claim ceiling:** "AOX-inhibitor-sensitive potency/respiration." Not AOX-gene
  causality (needs a knockout).
