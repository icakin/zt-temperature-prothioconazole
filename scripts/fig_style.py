"""fig_style.py — shared publication style (matplotlib side).
Mirror of fig_style.R — SAME hex values, SAME font, SAME sizes.
Validated with the dataviz palette validator: the 15/21/27 anchors pass all
categorical checks; full ramps are ORDERED scales (identity carried by order +
direct labels); grays are semantic neutrals (control), never a series color.
"""
import matplotlib as mpl

FONT = "DejaVu Sans"

# --- temperature: ordered cool->warm ramp (7 levels) --------------------------
TEMPS = [15, 18, 21, 24, 26, 27, 28]
TEMP_RAMP = ["#27519E", "#457BBE", "#64A5DE", "#C9A227",
             "#E08214", "#D64B21", "#A31621"]
TEMP_COL = dict(zip(TEMPS, TEMP_RAMP))
# RNA-seq uses three temperatures as identities — validated trio:
TEMP3 = {15: "#27519E", 21: "#64A5DE", 27: "#D64B21"}  # exact ramp members

# --- dose: control gray + single-hue purple ramp light->dark (ordered) --------
DOSES = ["Control", "0.06", "0.12", "0.25", "0.5", "1", "2", "4"]
DOSE_RAMP = ["#8C8C8C", "#A1DA38", "#4FC36B", "#21A585", "#25848E",
             "#33638D", "#433D84", "#471063"]
DOSE_COL = dict(zip(DOSES, DOSE_RAMP))

# --- control vs drug pair (bars, paired series) -------------------------------
CTRL = "#9A9A9A"
DRUG = "#433D84"          # ties to the dose ramp's dark end

# --- diverging (Bliss deviation, log2FC, NES): blue <-> red, neutral mid ------
DIV_CMAP = "RdBu_r"
DIV_LOW, DIV_MID, DIV_HIGH = "#2166AC", "#F7F7F7", "#B2182B"

# --- functional gene classes (RNA volcano highlights; fixed order) ------------
CLASS_COL = {
    "Efflux (ABC transporters)": "#7B3FA0",
    "Ergosterol / sterol":       "#1B7837",
    "Cytochrome-P450 detox":     "#D6604D",
    "Respiration (OXPHOS / AOX)":"#E08214",
    "Translation (ribosome)":    "#2F6BB3",
    "TCA cycle":                 "#8C510A",
}

GRAY_PT = "#C8C8C8"       # background points
INK = "#1A1A1A"           # primary text
INK2 = "#555555"          # secondary text

# --- sizes (figure built at 183 mm double-column width, 600 dpi export) -------
MM = 1 / 25.4
W2COL = 183 * MM          # 7.20 in
BASE = 7.0                # base font pt
LABEL = 10.0              # panel-label pt (bold)

def apply():
    mpl.rcParams.update({
        "font.family": FONT,
        "font.size": BASE,
        "axes.titlesize": BASE + 1, "axes.titleweight": "bold",
        "axes.labelsize": BASE, "axes.linewidth": 0.6,
        "axes.edgecolor": INK2, "axes.labelcolor": INK,
        "xtick.labelsize": BASE - 0.5, "ytick.labelsize": BASE - 0.5,
        "xtick.color": INK2, "ytick.color": INK2,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "legend.fontsize": BASE - 0.5, "legend.frameon": False,
        "axes.spines.top": False, "axes.spines.right": False,
        "figure.dpi": 110, "savefig.dpi": 600,
        "savefig.facecolor": "white", "figure.facecolor": "white",
        "pdf.fonttype": 42, "ps.fonttype": 42,
    })

def panel_label(ax, s, dx=-0.10, dy=1.06):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=LABEL,
            fontweight="bold", va="top", ha="left", color=INK)
