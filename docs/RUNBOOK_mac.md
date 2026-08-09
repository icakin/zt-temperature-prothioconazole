# zt_ipo323 — Mac runbook (the three remaining externals)

All three tasks below run on your Mac. Prereqs you already have from the
cauris work: the `temstapro_env_CPU` conda env, a TemStaPro clone, and Gab's
`etcGEMs` repo cloned with the package installed (`pip install -e .`).

Setup once: unzip `zt_ipo323_strain.zip` (from `X0123 copy/gem/`) into Gab's
repo so the folder sits at `etcGEMs/strains/zt_ipo323/`. Everything below
assumes you `cd` into the `etcGEMs` repo root, and that ZT points at the
strain:

    export ZT="strains/zt_ipo323"

---

## 1. TemStaPro — per-enzyme melting temperatures (~30–60 min CPU)

The enzyme sequences are already extracted: `$ZT/thermal/enzyme_sequences.faa`
(742 model enzymes). Mirroring your cauris `run_temstapro.sh`:

    conda activate temstapro_env_CPU
    TEMSTAPRO_DIR=/path/to/your/TemStaPro   # same clone you used for cauris
    mkdir -p "$ZT/thermal/temstapro_cache"

    python "$TEMSTAPRO_DIR/temstapro" \
        -f "$ZT/thermal/enzyme_sequences.faa" \
        -d "$TEMSTAPRO_DIR/ProtTrans" \
        -e "$ZT/thermal/temstapro_cache" \
        -t "$TEMSTAPRO_DIR" \
        --more-thresholds \
        --mean-output "$ZT/thermal/enzyme_tm_temstapro_raw.tsv"

    python "$ZT/scripts/temstapro_to_bestparams_zt.py" \
        --tsv "$ZT/thermal/enzyme_tm_temstapro_raw.tsv" \
        --faa "$ZT/thermal/enzyme_sequences.faa" \
        --out "$ZT/thermal/BestParamsTopt.csv"

(Unlike cauris there is no Z. tritici homology meltome, so the sequence
pseudo-Tm is used directly; the script flags this in TmTag. --center defaults
to 24 C, the measured growth optimum.)

Then edit `$ZT/strain.yaml`:
  - uncomment:  `enzyme_params: thermal/BestParamsTopt.csv`
  - change:     `thermal_model: mmrt`  ->  `thermal_model: unfolding`

And re-run the emergent curve:

    etcgem build --strain zt_ipo323
    etcgem tpc   --strain zt_ipo323

What to look for: the falling limb should sharpen and CTmax should move off
the grid edge (35 C) to something realistic (~29–31 C, near the measured
collapse). Compare against `$ZT/thermal/ExpGrowth.csv`.

## 2. DLTKcat — temperature-dependent kcats

Gab's pipeline generates the input for you:

    etcgem dltkcat prep --strain zt_ipo323 --tmin 5 --tmax 35 --n 11

That writes `$ZT/dltkcat/input.csv` (rxn_id, enzyme, substrate, SMILES,
sequence, temperatures). Run your DLTKcat clone on it (same as the eciML1515
workflow in Gab's docs/RUNBOOK.md), save predictions as
`$ZT/dltkcat/output.csv`, then:

    etcgem dltkcat parse --strain zt_ipo323 --temp-col Temp_C --kcat-col pred_log10kcat
    etcgem tpc --strain zt_ipo323

What to look for: rmax should move toward the measured 0.052 /h (currently
~4x high with uniform placeholder kcats).

## 3. CarveFungi check (~10 min, mostly download)

    curl -L -o /tmp/CarveFungi.tar.gz "https://zenodo.org/records/7413265/files/CarveFungi.tar.gz?download=1"
    tar -tzf /tmp/CarveFungi.tar.gz | grep -i -E "zymoseptoria|mycosphaerella|graminicola|tritici"

Empty output = no competing Z. tritici model in the 834-fungus collection,
and the novelty claim for ztGEM stands. If it DOES find one, extract it and
send it over — comparing it against ztGEM would take one session.

---

After 1 + 2, the strain is fully Layer-1 grounded and ready for Gab's
calibration stage (`etcgem calibrate`, his call on settings). Ship him the
folder plus `gem/ztGEM_v04_report.md` as the scientific context.
