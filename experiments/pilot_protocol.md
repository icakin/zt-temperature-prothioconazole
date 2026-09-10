# AOX-inhibitor pilot — full protocol

**Goal:** fix the nPG working dose, prove AOX blockade (azide-resistant respiration
collapses under nPG and SHAM), and measure autoxidation blanks — on 27 °C + 2 mg/L
prothioconazole cells. One SensorDish plate (24 vials).

**Solvent rule (Route A):** the drug stays in **0.04% DMSO** exactly as the main
study; the inhibitors go in **ethanol** (0.2% final, matched in every vial). DMSO
therefore never exceeds your 0.04%.

---

## 1. Reagent stocks

| Stock | Recipe | Solvent | Store |
|---|---|---|---|
| Sodium azide 1 M | 65 mg in 1 mL water | water | 4 °C, "toxic" |
| Prothioconazole 5 mg/mL | 5 mg in 1 mL DMSO | DMSO | −20 °C |
| nPG master 100 mM | 21 mg n-propyl gallate in 1 mL **ethanol** | ethanol | fresh, foil |
| SHAM master 500 mM | 76.5 mg N,2-dihydroxybenzamide in 1 mL **ethanol** | ethanol | fresh, foil |

> Remake the nPG/SHAM masters in **ethanol** (the DMSO versions were made before we
> fixed the 0.04% DMSO limit — keep them only as spares). nPG and SHAM are
> light-sensitive: foil-wrap, make on the day.

## 2. Inhibitor working series (ethanol; each is a 500× stock → 10 µL into 5 mL)

**nPG** (from the 100 mM ethanol master):

| stock | recipe | → final in vial |
|---|---|---|
| 12.5 mM | 12.5 µL master + 87.5 µL EtOH | 0.025 mM |
| 25 mM | 25 µL master + 75 µL EtOH | 0.05 mM |
| 50 mM | 50 µL master + 50 µL EtOH | 0.1 mM |
| 75 mM | 75 µL master + 25 µL EtOH | 0.15 mM |
| 100 mM | master, neat | 0.2 mM |

**SHAM** (from the 500 mM ethanol master):

| stock | recipe | → final in vial |
|---|---|---|
| 125 mM | 25 µL master + 75 µL EtOH | 0.25 mM |
| 250 mM | 50 µL master + 50 µL EtOH | 0.5 mM |
| 500 mM | master, neat | 1 mM |

## 3. Preculture (start ~4 days ahead — the real gate)

Inoculate from frozen stock into liquid YMS; grow **4 days, 21 °C, 200 rpm**; pass
through a 100 µm cell strainer to break clumps. Measure OD₆₀₀ on the day.

## 4. On the day — two master cups, then distribute

**Cup A — inoculated master (for the 13 cell vials), make 70 mL:**
- 70 mL YMS
- + 28 µL of the 5 mg/mL prothioconazole stock → 2 mg/L, 0.04% DMSO
- + preculture to OD₆₀₀ 0.0005 (35 µL of an OD 1.0 preculture; scale to your OD:
  µL = 0.0005 × 70 000 ÷ your OD)
- mix; distribute **4.985 mL per vial** to A1–A6, C1–C3, D1–D4

**Cup B — sterile master (for the 11 blank vials), make 60 mL:**
- 60 mL sterile YMS
- + 24 µL of the 5 mg/mL prothioconazole stock → 2 mg/L, 0.04% DMSO
- **no cells**
- mix; distribute **4.985 mL per vial** to B1–B6, C4–C6, D5–D6

## 5. One per-vial step — add 10 µL inhibitor (ethanol stock or plain ethanol)

- plain ethanol → A1, B1, D1, D2, D5
- nPG 12.5 / 25 / 50 / 75 / 100 mM → A2–A6 and B2–B6 (matching ladders)
- nPG 100 mM → D3, D6
- SHAM 125 / 250 / 500 mM → C1–C3 and C4–C6 (matching ladders)
- SHAM 500 mM → D4

## 6. Run

Seal → SensorDish baseline at 27 °C → inject **2.5 µL of 1 M azide** into D2–D6
through the septum → keep recording. Analyse the **initial high-O₂ window**, not the
whole depletion trace.

## 7. Full vial map

| Well | Bulk | Inhibitor (10 µL) | Final | Azide |
|---|---|---|---|---|
| A1 | inoculated | plain EtOH | nPG 0 | — |
| A2 | inoculated | 12.5 mM nPG | 0.025 mM | — |
| A3 | inoculated | 25 mM nPG | 0.05 mM | — |
| A4 | inoculated | 50 mM nPG | 0.1 mM | — |
| A5 | inoculated | 75 mM nPG | 0.15 mM | — |
| A6 | inoculated | 100 mM nPG | 0.2 mM | — |
| B1 | sterile | plain EtOH | nPG 0 | — |
| B2 | sterile | 12.5 mM nPG | 0.025 mM | — |
| B3 | sterile | 25 mM nPG | 0.05 mM | — |
| B4 | sterile | 50 mM nPG | 0.1 mM | — |
| B5 | sterile | 75 mM nPG | 0.15 mM | — |
| B6 | sterile | 100 mM nPG | 0.2 mM | — |
| C1 | inoculated | 125 mM SHAM | 0.25 mM | — |
| C2 | inoculated | 250 mM SHAM | 0.5 mM | — |
| C3 | inoculated | 500 mM SHAM | 1 mM | — |
| C4 | sterile | 125 mM SHAM | 0.25 mM | — |
| C5 | sterile | 250 mM SHAM | 0.5 mM | — |
| C6 | sterile | 500 mM SHAM | 1 mM | — |
| D1 | inoculated | plain EtOH | reference | none |
| D2 | inoculated | plain EtOH | — | 2.5 µL |
| D3 | inoculated | 100 mM nPG | 0.2 mM | 2.5 µL |
| D4 | inoculated | 500 mM SHAM | 1 mM | 2.5 µL |
| D5 | sterile | plain EtOH | — | 2.5 µL |
| D6 | sterile | 100 mM nPG | 0.2 mM | 2.5 µL |

## 8. Read-out

- **Dose-finding (A1–A6):** lowest nPG dose that blocks AOX (from D-arm) without
  denting growth vs A1, and with a clean matching blank (B-row).
- **Blockade proof (D-arm):** D2 shows azide-resistant (AOX) respiration; D3 (nPG)
  and D4 (SHAM) should both collapse it toward the D5/D6 abiotic floor.
- **Autoxidation:** subtract each cell trace's matched blank by fitting a background
  term; reject any dose whose blank O₂ loss rivals the biological signal.
