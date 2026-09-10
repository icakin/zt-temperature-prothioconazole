# AOX-inhibitor potency-restoration — main experiment protocol

**Question:** does blocking AOX (nPG, confirmed by SHAM) restore prothioconazole
potency at 27 °C but not at 15 °C? Headline signature = EC₅₀ falls with inhibitor
at 27 °C, barely moves at 15 °C (an inhibitor × temperature interaction).

**Design (one setup day):** 3 biological replicates = 3 independent precultures.
Each preculture → **one inoculated master → two plates, 15 °C + 27 °C** (same
culture at both temperatures — the pairing is the experiment). Each plate carries
the full **6-dose × 3-arm** design. Plus one **azide-validation** plate on dense
cells. **7 plates total.**

**Solvent (Route A):** drug in **0.04% DMSO** (2 µL of a dose stock per 5 mL);
inhibitors in **ethanol** (10 µL, 0.2% final, matched in every vial). DMSO never
exceeds 0.04%.

---

## 1. Reagent stocks

| Stock | Recipe | Solvent | Notes |
|---|---|---|---|
| Sodium azide 1 M | 65 mg / 1 mL water | water | 4 °C, "toxic" |
| nPG working (ethanol) | at pilot dose (see §2) | ethanol | foil, fresh |
| SHAM working (ethanol) | at pilot dose (see §2) | ethanol | foil, fresh |

**Prothioconazole dose ladder** — each dose delivered as **2 µL into 5 mL**
(= 0.04% DMSO). Make/confirm these DMSO stocks:

| Final in vial (mg/L) | DMSO stock (add 2 µL) |
|---|---|
| 0 | pure DMSO |
| 0.25 | 0.625 mg/mL |
| 0.5 | 1.25 mg/mL |
| 1 | 2.5 mg/mL |
| 2 | 5 mg/mL |
| 4 | 10 mg/mL |

## 2. Inhibitor working stocks (ethanol; single working dose each, set by the pilot)

Placeholders until the pilot fixes them — nPG ~0.1 mM, SHAM ~0.5 mM. Each is a
**500× ethanol stock** (10 µL → 5 mL):

- nPG: 0.1 mM final → **50 mM** ethanol stock (adjust to the pilot dose)
- SHAM: 0.5 mM final → **250 mM** ethanol stock (adjust to the pilot dose)
- "Control" arm: **plain ethanol**, 10 µL

## 3. Cultures (start ahead)

- **3 precultures** (biological replicates): frozen stock → YMS, 4 days, 21 °C,
  200 rpm, 100 µm strain. Measure OD₆₀₀ on the day.
- **1 dense flask for azide:** a shake-flask grown to **mid-exponential** at
  27 °C + 2 mg/L prothioconazole (dense enough that O₂ falls in ~1–2 h). This is
  NOT at OD 0.0005 — the acute assay needs cells already respiring hard.

## 4. On the day — per replicate: one master → two plates

Repeat for each of the 3 precultures.

**Inoculated master, make 190 mL:**
- 190 mL YMS
- + preculture to OD₆₀₀ 0.0005 (95 µL of an OD 1.0 preculture; scale:
  µL = 0.0005 × 190 000 ÷ your OD)
- **no drug, no inhibitor** (both vary per vial)
- mix; use it to fill the **18 culture wells of the 15 °C plate and the 18 of the
  27 °C plate** (4.99 mL/well)

**Sterile master, make 65 mL:**
- 65 mL YMS, no cells
- fill the **6 blank wells of each plate** (4.99 mL/well)

## 5. Per-vial additions (both plates identical)

Every vial: **4.99 mL master + 2 µL drug (that dose's stock) + 10 µL inhibitor**
(plain ethanol / nPG / SHAM). Culture wells use inoculated master; blanks use
sterile master.

## 6. Plate map — same for the 15 °C and 27 °C plates

Rows = inhibitor arm; columns = prothioconazole dose.

| | 1 (0) | 2 (0.25) | 3 (0.5) | 4 (1) | 5 (2) | 6 (4 mg/L) |
|---|---|---|---|---|---|---|
| **A control** | ethanol | ethanol | ethanol | ethanol | ethanol | ethanol |
| **B nPG** | +nPG | +nPG | +nPG | +nPG | +nPG | +nPG |
| **C SHAM** | +SHAM | +SHAM | +SHAM | +SHAM | +SHAM | +SHAM |
| **D blanks** | ctrl 0 | ctrl 4 | nPG 0 | nPG 4 | SHAM 0 | SHAM 4 |

Rows A–C are cells (inoculated master); each cell = 2 µL of the column's drug
stock + 10 µL of the row's inhibitor. Row D = cell-free blanks (sterile master),
same chemistry, no cells — the autoxidation background for the 0 and 4 mg/L ends
of each arm.

## 7. Azide-validation plate (dense cells, acute)

Fill from the **dense 27 °C + drug flask** (not OD 0.0005). 6 vials:

| Well | Cells | Inhibitor (10 µL) | Inject after baseline |
|---|---|---|---|
| A1 | dense | plain ethanol | nothing (total respiration) |
| A2 | dense | plain ethanol | 2.5 µL 1 M azide → azide-resistant = AOX |
| A3 | dense | nPG | 2.5 µL azide → should collapse A2's residual |
| A4 | dense | SHAM | 2.5 µL azide → should also collapse it |
| A5 | sterile | plain ethanol | 2.5 µL azide → abiotic blank |
| A6 | sterile | nPG | 2.5 µL azide → nPG blank |

All at 2 mg/L drug (2 µL drug stock each). Inject azide once the O₂ trace is in
**steep decline** (dense cells respiring hard), read the slope drop.

## 8. Run & handling

- Seal all plates; 15 °C and 27 °C plates to their readers; run **~4 days** to full
  O₂ depletion for the growth curves.
- **Keep everything dark** (foil / dark incubator): nPG and SHAM photodegrade, and
  over a 4-day chronic run the effective dose drifts down — the blanks and the
  azide plate are your checks that inhibition held.
- Azide plate: acute, ~1–2 h, done whenever the dense flask is ready.

## 9. Read-out / analysis

- Fit each DO trace (manuscript oxygen model) → growth rate r + biomass-specific
  respiration; correct each with its matched blank (fit a background term, don't
  subtract one slope; reject a dose whose blank rivals the signal).
- **Primary:** Hill fit per (temperature × arm) → EC₅₀. Signature = EC₅₀(27 °C,
  +nPG) ≪ EC₅₀(27 °C, control), EC₅₀(15 °C) little changed; report the
  inhibitor × temperature interaction on log EC₅₀. SHAM must reproduce the nPG
  rescue.
- **Secondary:** is the drug-induced extra respiration nPG/SHAM-sensitive?
  (compare fitted respiration ± inhibitor per temperature.)
- **Azide plate:** A2 shows AOX-carried respiration; A3/A4 should collapse it —
  the blockade proof.
- Replication: culture/day = biological replicate (n = 3 here); mixed model with a
  random intercept per preculture, or paired per-replicate contrasts. This is a
  first, effect-size block — extend with more precultures on later days to power
  the EC₅₀ shift fully.
