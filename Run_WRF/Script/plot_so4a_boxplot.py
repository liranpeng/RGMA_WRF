"""
plot_so4a_boxplot.py

Estimate so4aj (acc) and so4ai (ait) mass mixing ratios [kg/kg] from the
MERRA-2 Nt (>100 nm) median values read from the box-plot figure, then
reproduce a similar box plot.

Formulas
--------
  so4aj = 1e6 * NTOT [cm⁻³] * rho_aer * (pi/6) * Dacc^3
          * exp(9/2 * ln(sigma_acc)^2) / rho_air

  so4ai = 2 * 1e6 * NTOT [cm⁻³] * rho_aer * (pi/6) * Dait^3
          * exp(9/2 * ln(sigma_nuc)^2) / rho_air

Both are linear in NTOT, so the full box statistics transform by the same
constant factor.  Box stats (whisker_lo, Q1, median, Q3, whisker_hi) for
the MERRA-2 red boxes are estimated visually from the original log-scale
figure.
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Physical constants ────────────────────────────────────────────────────────
pi        = np.pi
rho_aer   = 1.8e3    # kg m⁻³  (SO4 aerosol density)
rho_air   = 1.225    # kg m⁻³  (standard air density)

# Mode parameters (SORGAM defaults)
Dacc      = 0.16e-6  # m   accumulation mode mean diameter
Dait      = 0.04e-6  # m   Aitken mode mean diameter
sigma_acc = 2.00     # geometric standard deviation, acc
sigma_nuc = 1.70     # geometric standard deviation, ait

# Pre-compute log-normal moment factor  exp(9/2 * ln(σ)²)
es36_acc = np.exp(4.5 * np.log(sigma_acc) ** 2)   # ≈ 8.689
es36_nuc = np.exp(4.5 * np.log(sigma_nuc) ** 2)   # ≈ 3.550

# Conversion factors [kg kg⁻¹ per cm⁻³]
#   so4aj = NTOT_cm3 * fac_acc
#   so4ai = NTOT_cm3 * fac_nuc  (the ×2 is already folded in)
fac_acc = 1e6 * rho_aer * (pi / 6) * Dacc**3 * es36_acc / rho_air
fac_nuc = 2 * 1e6 * rho_aer * (pi / 6) * Dait**3 * es36_nuc / rho_air

print(f"es36_acc = {es36_acc:.4f}")
print(f"es36_nuc = {es36_nuc:.4f}")
print(f"fac_acc  = {fac_acc:.4e} kg/kg per cm⁻³")
print(f"fac_nuc  = {fac_nuc:.4e} kg/kg per cm⁻³")
print(f"so4aj / so4ai ratio ≈ {fac_acc / fac_nuc:.1f}")

# ── MERRA-2 Nt(>100 nm) box stats read from figure [cm⁻³] ────────────────────
# Columns: [whisker_lo, Q1, median, Q3, whisker_hi]
# Medians are the labeled values; quartiles / whiskers are estimated visually
# from the log-scale red boxes in the original figure.

dates = [
    "2017-07-11", "2017-07-12", "2017-07-13", "2017-07-15",
    "2017-07-17", "2017-07-18", "2017-07-19", "2017-07-20",
]

nt_boxes = np.array([
    [ 85,  100, 117.2, 130, 152],   # 2017-07-11  median labeled 117.2
    [132,  172, 210.6, 252, 315],   # 2017-07-12  median labeled 210.6
    [128,  158, 198.8, 238, 292],   # 2017-07-13  median labeled 198.8
    [ 47,   53,  62.2,  72,  88],   # 2017-07-15  median labeled  62.2
    [ 97,  112, 128.4, 146, 168],   # 2017-07-17  median labeled 128.4
    [ 90,  103, 116.6, 130, 154],   # 2017-07-18  median labeled 116.6
    [ 53,   63,  73.9,  84,  97],   # 2017-07-19  median labeled  73.9
    [ 97,  113, 130.7, 150, 175],   # 2017-07-20  median labeled 130.7
])

# Apply linear conversion
so4aj_boxes = nt_boxes * fac_acc   # kg/kg
so4ai_boxes = nt_boxes * fac_nuc   # kg/kg

# ── Plotting ──────────────────────────────────────────────────────────────────
fig, (ax, ax2) = plt.subplots(
    2, 1, figsize=(16, 8),
    gridspec_kw={"height_ratios": [3, 1]},
)

n       = len(dates)
x       = np.arange(n)
bw      = 0.28                          # box full-width
offsets = [-0.32, 0.32]                 # so4aj left, so4ai right
colors  = ["#1565C0", "#C62828"]        # blue=acc, red=ait
labels  = ["so4aj (acc mode)", "so4ai (ait mode)"]
all_boxes = [so4aj_boxes, so4ai_boxes]

for boxes, offset, color, label in zip(all_boxes, offsets, colors, labels):
    for i in range(n):
        wl, q1, med, q3, wh = boxes[i]
        cx = x[i] + offset             # center x of this box

        # IQR rectangle
        rect = mpatches.FancyBboxPatch(
            (cx - bw / 2, q1), bw, q3 - q1,
            boxstyle="square,pad=0",
            facecolor=color, alpha=0.65,
            edgecolor="black", linewidth=0.9,
        )
        ax.add_patch(rect)

        # Median line
        ax.plot([cx - bw / 2, cx + bw / 2], [med, med],
                color="black", linewidth=1.8, solid_capstyle="butt")

        # Whiskers
        ax.plot([cx, cx], [wl, q1],  color="black", linewidth=0.9)
        ax.plot([cx, cx], [q3, wh],  color="black", linewidth=0.9)
        cap = bw * 0.35
        ax.plot([cx - cap, cx + cap], [wl, wl], color="black", linewidth=0.9)
        ax.plot([cx - cap, cx + cap], [wh, wh], color="black", linewidth=0.9)

        # Median value label
        ax.text(cx, wh * 1.12, f"{med:.2e}",
                ha="center", va="bottom", fontsize=6.5,
                color=color, fontweight="bold")

# Axes formatting
ax.set_yscale("log")
ax.set_xlim(-0.7, n - 0.3)
ax.set_xticks(x)
ax.set_xticklabels(dates, rotation=30, ha="right", fontsize=10)
ax.set_ylabel("Sulfate mass mixing ratio [kg kg⁻¹]", fontsize=11)
ax.set_title(
    "Box summary: 2017-07-11 to 2017-07-20\n"
    "so4aj (acc) and so4ai (ait) estimated from MERRA-2 Nt(>100 nm)",
    fontsize=12,
)
ax.grid(True, axis="y", which="both", alpha=0.3, linestyle="--")

# Legend
legend_handles = [
    mpatches.Patch(facecolor=colors[0], alpha=0.65, edgecolor="black",
                   label=f"so4aj (acc)  fac={fac_acc:.3e} kg/kg per cm⁻³"),
    mpatches.Patch(facecolor=colors[1], alpha=0.65, edgecolor="black",
                   label=f"so4ai (ait)  fac={fac_nuc:.3e} kg/kg per cm⁻³"),
]
ax.legend(handles=legend_handles, loc="upper right", fontsize=9)

# Annotation box with formula
formula_text = (
    r"$\mathrm{so4aj} = 10^6 \cdot N_t \cdot \rho_{aer} \cdot \frac{\pi}{6} \cdot D_{acc}^3 "
    r"\cdot e^{\frac{9}{2}\ln^2\sigma_{acc}} / \rho_{air}$"
    "\n"
    r"$\mathrm{so4ai} = 2\times10^6 \cdot N_t \cdot \rho_{aer} \cdot \frac{\pi}{6} \cdot D_{ait}^3 "
    r"\cdot e^{\frac{9}{2}\ln^2\sigma_{ait}} / \rho_{air}$"
    "\n"
    rf"$\rho_{{aer}}={rho_aer:.0f}$ kg m⁻³,  $D_{{acc}}={Dacc*1e6:.2f}$ μm,"
    rf"  $\sigma_{{acc}}={sigma_acc}$,  $D_{{ait}}={Dait*1e6:.2f}$ μm,  $\sigma_{{ait}}={sigma_nuc}$"
)
ax.text(0.01, 0.03, formula_text, transform=ax.transAxes,
        fontsize=7.5, va="bottom", ha="left",
        bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow",
                  edgecolor="gray", alpha=0.85))

# ── Bottom panel: stacked bar – mass fraction between so4aj and so4ai ─────────
saj_med = so4aj_boxes[:, 2]          # median so4aj per date
sai_med = so4ai_boxes[:, 2]          # median so4ai per date
total   = saj_med + sai_med

frac_aj = saj_med / total            # acc fraction
frac_ai = sai_med / total            # ait fraction

bar_kw = dict(width=0.55, edgecolor="black", linewidth=0.7)
ax2.bar(x, frac_aj, color="#1565C0", alpha=0.65, label="so4aj (acc)", **bar_kw)
ax2.bar(x, frac_ai, bottom=frac_aj,  color="#C62828", alpha=0.65, label="so4ai (ait)", **bar_kw)

ax2.set_xlim(-0.7, n - 0.3)
ax2.set_xticks(x)
ax2.set_xticklabels(dates, rotation=30, ha="right", fontsize=10)
ax2.set_yticks([0, 0.5, 1.0])
ax2.set_yticklabels(["0", "0.5", "1"], fontsize=9)
ax2.set_ylabel("Mass\nfraction", fontsize=9)
ax2.set_ylim(0, 1)
ax2.grid(True, axis="y", alpha=0.3, linestyle="--")
ax2.legend(loc="upper right", fontsize=8, ncol=2)

plt.tight_layout()
out_path = "Run_WRF/Script/so4a_boxplot.png"
plt.savefig(out_path, dpi=150, bbox_inches="tight")
print(f"Figure saved → {out_path}")
plt.show()

# ── Print summary table ───────────────────────────────────────────────────────
print(f"\n{'Date':<13} {'Nt_med':>8} {'so4aj_med':>14} {'so4ai_med':>14}"
      f"  {'so4aj range':>22}  {'so4ai range':>22}")
print("─" * 100)
for i, date in enumerate(dates):
    nt  = nt_boxes[i]
    saj = so4aj_boxes[i]
    sai = so4ai_boxes[i]
    print(f"{date:<13} {nt[2]:>8.1f} {saj[2]:>14.3e} {sai[2]:>14.3e}"
          f"  [{saj[0]:.3e} – {saj[4]:.3e}]  [{sai[0]:.3e} – {sai[4]:.3e}]")
