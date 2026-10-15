import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
from matplotlib.path import Path
from matplotlib.patches import PathPatch

plt.rcParams.update({
    "font.size": 15, "axes.labelsize": 15, "axes.titlesize": 15,
    "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 15,
    "mathtext.fontset": "cm", "font.family": "serif", "font.serif": ["cmr10"],
    "axes.formatter.use_mathtext": True, "axes.unicode_minus": False,
    "svg.fonttype": "path",
})

FLIP, LOC = "#fcd5a8", "#c4e3e1"
UP, DN, RESP = "#b01c1c", "#1c3ab0", "#9b1c1c"
FS = 15

def vbrace(ax, x, y0, y1, w=0.08, color="k", lw=1.2, side=1):
    """curly brace along y from y0 to y1 at x, opening toward +x if side=1"""
    ym = 0.5*(y0+y1); s = side
    verts = [(x, y0), (x+s*w, y0), (x+s*w, ym), (x+2*s*w, ym),
             (x+s*w, ym), (x+s*w, y1), (x, y1)]
    codes = [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.LINETO,
             Path.MOVETO, Path.CURVE3, Path.CURVE3]
    # build as two quadratic halves
    p = Path([(x, y0), (x+s*w, y0), (x+s*w, y0+0.35*(ym-y0)),
              (x+s*w, ym-0.12*(ym-y0)), (x+2*s*w, ym),
              (x+s*w, ym+0.12*(y1-ym)), (x+s*w, ym+0.65*(y1-ym)),
              (x+s*w, y1), (x, y1)],
             [Path.MOVETO, Path.CURVE3, Path.CURVE3,
              Path.LINETO, Path.LINETO,
              Path.LINETO, Path.LINETO,
              Path.CURVE3, Path.CURVE3])
    ax.add_patch(PathPatch(p, fc="none", ec=color, lw=lw, capstyle="round", joinstyle="round"))

# ======================= Figure (a): 6 x 4 in =======================
pattern = np.array([1, -1, -1, 1, -1])
fig = plt.figure(figsize=(6, 4))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 6); ax.set_ylim(0, 4); ax.axis("off")

t0, T, T1 = 0.40, 1.55, 0.5          # y of t=0, period, T1 (inches, time runs upward)
yT = lambda n: t0 + n*T
xa, xb = 0.62, 1.10                  # drive block edges / time axis at xa
xc = 2.55                            # chain centre


for n in range(2):
    ya, yb, yc = yT(n), yT(n)+T1, yT(n+1)
    ax.add_patch(Rectangle((xa, ya), xb-xa, yb-ya, fc=FLIP, ec="0.35", lw=0.9))
    ax.add_patch(Rectangle((xa, yb), xb-xa, yc-yb, fc=LOC, ec="0.35", lw=0.9))
    ax.text((xa+xb)/2, (ya+yb)/2, r"$\hat H_1$", ha="center", va="center", fontsize=FS)
    ax.text((xa+xb)/2, (yb+yc)/2, r"$\hat H_2$", ha="center", va="center", fontsize=FS)

ax.annotate("", xy=(xa, 3.9), xytext=(xa, t0-0.22),
            arrowprops=dict(arrowstyle="-|>", lw=1.4, color="k", mutation_scale=14))
ax.text(xa+0.1, 3.86, r"$t$", fontsize=FS, va="center")
for n, lab in enumerate([r"$0$", r"$T$", r"$2T$"]):
    ax.plot([xa-0.07, xa+0.05], [yT(n)]*2, color="k", lw=1.4)
    ax.text(xa-0.11, yT(n), lab, ha="right", va="center", fontsize=FS)

vbrace(ax, xb+0.05, yT(0), yT(0)+T1, w=0.06)
ax.text(xb+0.22, yT(0)+T1/2, r"$T_1$", va="center", fontsize=FS)
vbrace(ax, xb+0.05, yT(0)+T1, yT(1), w=0.06)
ax.text(xb+0.22, yT(0)+T1+(T-T1)/2, r"$T_2$", va="center", fontsize=FS)

for n in range(3):
    y = yT(n); sgn = 1 if n % 2 == 0 else -1
    x0 = xb+0.62 if n == 0 else xb
    ax.plot([x0, xc-0.62], [y, y], ls=(0, (1, 2)), color="0.5", lw=1)
    ax.plot([xc-0.6, xc+0.6], [y, y], color="0.65", lw=1)
    for k, s in enumerate(sgn*pattern):
        xx = xc - 0.48 + 0.24*k
        ax.plot(xx, y, "o", ms=3, color="0.5")
        col = UP if s > 0 else DN
        ax.annotate("", xy=(xx, y+0.24*s), xytext=(xx, y-0.2*s),
                    arrowprops=dict(arrowstyle="-|>", lw=2.0, color=col,
                                    mutation_scale=11, shrinkA=0, shrinkB=0))
    lab = r"$|\psi_0\rangle$" if sgn > 0 else r"$\hat X|\psi_0\rangle$"
    ax.text(xc+0.68, y, lab, va="center", fontsize=FS)

for n in range(2):
    ax.annotate("", xy=(xc, yT(n+1)-0.33), xytext=(xc, yT(n)+0.33),
                arrowprops=dict(arrowstyle="-|>", lw=1.6, color="0.3", mutation_scale=14))
    ax.text(xc+0.08, (yT(n)+yT(n+1))/2, r"$\hat U_F$", va="center", fontsize=FS, color="0.2")

# period-2T brace
vbrace(ax, 3.92, yT(0)-0.05, yT(2)+0.05, w=0.09, color=RESP, lw=1.6)
ax.text(4.2, (yT(0)+yT(2))/2, r"period $2T$", rotation=-90, ha="center", va="center",
        fontsize=FS, color=RESP)

# definitions (right column)
xd = 4.45
ax.add_patch(Rectangle((xd, 3.42), 0.22, 0.2, fc=FLIP, ec="0.35", lw=0.9))
ax.text(xd+0.3, 3.52, "spin flip", va="center", fontsize=FS)
ax.add_patch(Rectangle((xd, 3.02), 0.22, 0.2, fc=LOC, ec="0.35", lw=0.9))
ax.text(xd+0.3, 3.12, "localization", va="center", fontsize=FS)
ax.text(xd, 2.55, r"$gT_1=\frac{\pi}{2}(1-\epsilon)$", va="center", fontsize=FS)
ax.text(xd, 2.0, r"$\hat X=\prod_i\hat\sigma^x_i$", va="center", fontsize=FS)
ax.text(xd, 1.45, r"$\hat U_F=$", va="center", fontsize=FS)
ax.text(xd+0.12, 1.0, r"$e^{-i\hat H_2T_2}e^{-i\hat H_1T_1}$", va="center", fontsize=FS)
ax.text(xd, 0.45, r"$\omega_D=2\pi/T$", va="center", fontsize=FS)

for ext, kw in [("svg", {}), ("pdf", {}), ("jpg", dict(dpi=600, pil_kwargs={"quality": 95}))]:
    fig.savefig(f"dtc_fig_a.{ext}", **kw)
plt.close(fig)

# ======================= Figure (b)+(c): 7.5 x 6 in =======================
fig, (axb, axc) = plt.subplots(2, 1, figsize=(7.5, 6), constrained_layout=True)
n = np.arange(25); eps = 0.1
m_dtc = np.cos(np.pi*n)*(0.88 + 0.12*np.exp(-n/3))
m_free = np.cos(np.pi*n)*np.cos(np.pi*eps*n)
axb.axhline(0, color="0.6", ls="--", lw=0.8)
axb.plot(n, m_dtc, "-o", color=UP, ms=5, lw=1.1, label=r"DTC: $J_i\neq0$")
axb.plot(n, m_free, "-o", color=DN, ms=5, lw=1.1, mfc="white", alpha=0.85, label=r"free: $J_i=0$")
axb.set_xlabel(r"$n=t/T$"); axb.set_ylabel(r"$\langle\hat\sigma^z_i(nT)\rangle$")
axb.set_xlim(-0.5, 24.5); axb.set_ylim(-1.2, 1.2)
axb.set_xticks(range(0, 25, 4)); axb.set_yticks([-1, 0, 1])
axb.legend(loc="lower center", bbox_to_anchor=(0.5, 1.0), ncol=2, frameon=False,
           title=r"$\epsilon=0.1$", title_fontsize=15, borderaxespad=0.1)
axb.text(-0.13, 1.12, r"$\mathbf{(a)}$", transform=axb.transAxes, fontsize=FS)

w = np.linspace(0, 1.05, 3000); g = 0.009
lor = lambda w0: 1/(1+((w-w0)/g)**2)
axc.plot(w, lor(0.5), color=UP, lw=1.3, label="DTC")
axc.plot(w, 0.5*lor(0.45)+0.5*lor(0.55), color=DN, lw=1.3, alpha=0.85, label="free")
axc.axvline(1, color="0.5", ls="--", lw=1)
axc.text(0.985, 1.12, r"drive $\omega_D$", ha="right", va="center", color="0.4", fontsize=FS)
axc.annotate(r"peak at $\omega_D/2$", xy=(0.512, 0.88), xytext=(0.6, 0.88), va="center",
             color=UP, fontsize=FS, arrowprops=dict(arrowstyle="->", color=UP, lw=1.2))
axc.annotate("", xy=(0.44, 0.58), xytext=(0.56, 0.58),
             arrowprops=dict(arrowstyle="<->", color=DN, lw=1.2, shrinkA=0, shrinkB=0))
axc.text(0.42, 0.58, r"$\epsilon\,\omega_D$", ha="right", va="center", color=DN, fontsize=FS)
axc.set_xlim(0, 1.06); axc.set_ylim(0, 1.25)
axc.set_xticks([0, 0.25, 0.5, 0.75, 1])
axc.set_xticklabels([r"$0$", r"$\frac{1}{4}$", r"$\frac{1}{2}$", r"$\frac{3}{4}$", r"$1$"])
axc.set_yticks([])
axc.set_xlabel(r"$\omega/\omega_D$"); axc.set_ylabel(r"$\mathrm{FFT}[\hat\sigma^z_i(nT)]$")
axc.legend(loc="upper left", frameon=False)
axc.text(-0.13, 1.03, r"$\mathbf{(b)}$", transform=axc.transAxes, fontsize=FS)

for ext, kw in [("svg", {}), ("pdf", {}), ("jpg", dict(dpi=600, pil_kwargs={"quality": 95}))]:
    fig.savefig(f"dtc_fig_bc.{ext}", **kw)
plt.close(fig)
print("done")
