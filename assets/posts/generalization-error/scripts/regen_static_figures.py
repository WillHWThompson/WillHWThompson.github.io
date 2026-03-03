"""
Regenerate bias_variance_tradeoff.png and double_descent.png
using the unified muted-pastel palette from the interactive viz.

Palette (matches generative-model-viz.js / interactive-plots.js):
  ROSE    #c8a0a0   dusty rose   — high-bias / overfit signal
  SAGE    #7a9db5   steel sage   — variance / test error
  OCHRE   #c4974a   warm ochre   — threshold / teacher direction
  DARK    #3a3530   near-black   — total / primary line
  MIST    #b3c4cc   light steel  — secondary
  BG      #fafaf8   off-white    — figure background
  GRID    #e5dfd9   warm gray    — grid lines
  SPINE   #d0c8c0   light warm   — axes spines
  TEXT    #5a5550   warm gray    — labels, ticks
"""

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Palette ────────────────────────────────────────────────────────────────────
ROSE  = "#c8a0a0"
SAGE  = "#7a9db5"
OCHRE = "#c4974a"
DARK  = "#3a3530"
MIST  = "#b3c4cc"
BG    = "#fafaf8"
GRID  = "#e5dfd9"
SPINE = "#d0c8c0"
TEXT  = "#5a5550"
CREAM = "#ede8e3"

# ── rcParams ───────────────────────────────────────────────────────────────────
mpl.rcParams.update({
    "figure.facecolor":   BG,
    "axes.facecolor":     BG,
    "savefig.facecolor":  BG,
    "axes.edgecolor":     SPINE,
    "axes.linewidth":     0.75,
    "axes.spines.top":    False,
    "axes.spines.right":  False,
    "axes.grid":          True,
    "grid.color":         GRID,
    "grid.linewidth":     0.75,
    "grid.linestyle":     "-",
    "xtick.color":        TEXT,
    "ytick.color":        TEXT,
    "xtick.labelsize":    11,
    "ytick.labelsize":    11,
    "xtick.direction":    "out",
    "ytick.direction":    "out",
    "xtick.major.size":   4,
    "ytick.major.size":   4,
    "text.color":         DARK,
    "axes.labelcolor":    TEXT,
    "axes.labelsize":     12,
    "axes.titlesize":     13,
    "axes.titlepad":      10,
    "legend.frameon":     False,
    "legend.fontsize":    11,
    "font.family":        "sans-serif",
    "font.sans-serif":    ["Arial", "DejaVu Sans", "Helvetica Neue"],
    "figure.dpi":         150,
})

ARROW = dict(arrowstyle="->", lw=1.2, color=TEXT)


# ══════════════════════════════════════════════════════════════════════════════
# Figure 1 — Bias–Variance Tradeoff
# ══════════════════════════════════════════════════════════════════════════════

h        = np.linspace(0, 10, 400)
bias2    = 1.5 * np.exp(-0.4 * h)
variance = 0.05 * np.exp(0.5 * (h - 4))
total    = bias2 + variance
opt_idx  = np.argmin(total)
h_opt    = h[opt_idx]

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(h, bias2,    color=ROSE,  lw=2.5, label="Bias²")
ax.plot(h, variance, color=SAGE,  lw=2.5, label="Variance")
ax.plot(h, total,    color=DARK,  lw=2.0, ls="--", label="Total Error", alpha=0.8)

ax.axvline(h_opt, color=OCHRE, ls=":", lw=1.5)
ax.text(h_opt + 0.25, 2.1, "Optimal\nComplexity",
        color=OCHRE, fontsize=10, va="center", style="italic")

ax.annotate("Underfitting\n(high bias)",
            xy=(1.0, bias2[np.abs(h - 1.0).argmin()]),
            xytext=(2.2, 1.6),
            arrowprops={**ARROW, "color": ROSE},
            color=ROSE, fontsize=10, ha="center")

ax.annotate("Overfitting\n(high variance)",
            xy=(8.0, variance[np.abs(h - 8.0).argmin()]),
            xytext=(6.6, 2.5),
            arrowprops={**ARROW, "color": SAGE},
            color=SAGE, fontsize=10, ha="center")

ax.set_xlabel("Model Complexity")
ax.set_ylabel("Error")
ax.set_title("Bias–Variance Tradeoff")
ax.set_xlim(0, 10)
ax.set_ylim(0, 3.0)
ax.set_xticks([])          # complexity has no units
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig("figures/bias_variance_tradeoff.png", dpi=150)
plt.close(fig)
print("Saved figures/bias_variance_tradeoff.png")


# ══════════════════════════════════════════════════════════════════════════════
# Figure 2 — Double Descent
# ══════════════════════════════════════════════════════════════════════════════

h   = np.linspace(0.1, 10, 600)
thr = 4.0

bias2    = 1.5 * np.exp(-0.4 * h)
variance = 0.05 * np.exp(0.5 * (h - 4))
classical = bias2 + variance

spike      = 3.2 * np.exp(-((h - thr) / 0.09) ** 2)
over_error = 0.55 * np.exp(-0.65 * (h - thr)) + 0.18
test_error = np.where(h < thr, classical, over_error) + spike
train_error = 0.08 * np.exp(-0.35 * h)

fig, ax = plt.subplots(figsize=(8, 5))

ax.plot(h, test_error,  color=SAGE,  lw=2.5, ls="--", label="Test Error")
ax.plot(h, train_error, color=ROSE,  lw=2.5,           label="Train Error")

ax.axvline(thr, color=OCHRE, ls=":", lw=1.5)
ax.text(thr + 0.2, np.max(test_error) * 0.68, "Interpolation\nThreshold",
        color=OCHRE, fontsize=10, va="center", style="italic")

# Regime shading
ax.axvspan(0.1, thr, alpha=0.04, color=ROSE,  zorder=0)
ax.axvspan(thr, 10,  alpha=0.04, color=SAGE,  zorder=0)

y_top = np.max(test_error) * 1.18

ax.annotate("Classical\nregime",
            xy=(2.0, test_error[np.abs(h - 2.0).argmin()]),
            xytext=(1.2, y_top * 0.82),
            arrowprops={**ARROW, "color": ROSE},
            color=ROSE, fontsize=10, ha="center")

ax.annotate("Peak risk\n(p ≈ n)",
            xy=(thr, test_error[np.abs(h - thr).argmin()]),
            xytext=(thr + 1.1, y_top * 1.0),
            arrowprops={**ARROW, "color": OCHRE},
            color=OCHRE, fontsize=10, ha="left")

ax.annotate("Overparameterized\nregime",
            xy=(7.5, test_error[np.abs(h - 7.5).argmin()]),
            xytext=(6.8, y_top * 0.72),
            arrowprops={**ARROW, "color": SAGE},
            color=SAGE, fontsize=10, ha="center")

ax.set_xlabel("Model Complexity  (p / n →)")
ax.set_ylabel("Error")
ax.set_title("Double Descent")
ax.set_xlim(0.1, 10)
ax.set_ylim(0, y_top)
ax.set_xticks([])
ax.legend(loc="upper right")
fig.tight_layout()
fig.savefig("figures/double_descent.png", dpi=150)
plt.close(fig)
print("Saved figures/double_descent.png")
