import numpy as np
from scipy.linalg import expm
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------- Free-fermion (Majorana) description ----------
# U_F = exp(-i sum_j g_j sx_j) exp(-i sum_j J_j sz_j sz_{j+1}),  open chain of L spins
# Heisenberg: gamma -> R gamma with R = exp(A_K) exp(A_I)
def R_matrix(J, g):
    L = len(g)
    AK = np.zeros((2 * L, 2 * L)); AI = np.zeros((2 * L, 2 * L))
    for j in range(L):              # sx_j = i g_{2j-1} g_{2j}  (0-based: 2j, 2j+1)
        AK[2 * j, 2 * j + 1] = 2 * g[j]; AK[2 * j + 1, 2 * j] = -2 * g[j]
    for j in range(L - 1):          # sz_j sz_{j+1} = i g_{2j} g_{2j+1} (0-based: 2j+1, 2j+2)
        AI[2 * j + 1, 2 * j + 2] = 2 * J[j]; AI[2 * j + 2, 2 * j + 1] = -2 * J[j]
    return expm(AK) @ expm(AI)

def edge_autocorr(R, nmax):
    v = np.zeros(R.shape[0]); v[0] = 1.0
    out = []
    w = v.copy()
    for n in range(nmax + 1):
        out.append(w[0])           # (R^n)_{11}
        w = R.T @ w if False else w @ R  # row vector times R gives (R^n)_{1,:}
    return np.array(out)

def pi_bound(R):
    # weight of gamma_1 on the eigenvalue -1 eigenspace of R (Mazur-type bound)
    ev, V = np.linalg.eig(R)
    idx = np.where(np.abs(ev + 1) < 1e-6)[0]
    if len(idx) == 0:
        # near -1 (exponentially split pair in finite chain)
        idx = np.argsort(np.abs(ev + 1))[:2]
    Q, _ = np.linalg.qr(V[:, idx])
    return np.sum(np.abs(Q[0, :]) ** 2)

# ---------- Panel (a): single-particle quasienergies vs epsilon ----------
L = 40; J0 = np.pi / 8
eps_list = np.linspace(0, 0.6, 121)
qe = []
for eps in eps_list:
    R = R_matrix(J0 * np.ones(L - 1), (np.pi / 2) * (1 - eps) * np.ones(L))
    th = np.angle(np.linalg.eigvals(R))
    qe.append(np.abs(th))
qe = np.array(qe)

# ---------- Panel (b): edge autocorrelation and the bound ----------
L = 150; nmax = 200
edge_data = {}
for eps in [0.0, 0.1, 0.2, 0.3]:
    R = R_matrix(J0 * np.ones(L - 1), (np.pi / 2) * (1 - eps) * np.ones(L))
    C = edge_autocorr(R, nmax)
    edge_data[eps] = (C, pi_bound(R))
    sgn = (-1) ** np.arange(nmax + 1)
    print(f"eps={eps}: bound={edge_data[eps][1]:.4f}, late-time avg={np.mean((C*sgn)[100:]):.4f}")

# ---------- Panel (c): bulk spin autocorrelation (determinant formula) ----------
def bulk_corr(R, j, nmax):
    # infinite-temperature <sz_j(nT) sz_j> = det[(R^n)_{S,S}], S = first 2j-1 Majoranas
    S = 2 * j - 1; out = []; Rn = np.eye(R.shape[0])
    for n in range(nmax + 1):
        out.append(np.linalg.det(Rn[:S, :S])); Rn = Rn @ R
    return np.array(out)

Lc = 80; jc = 40; nmax_c = 200; eps_c = 0.1
sgn_c = (-1) ** np.arange(nmax_c + 1)
Rc = R_matrix(J0 * np.ones(Lc - 1), (np.pi / 2) * (1 - eps_c) * np.ones(Lc))
clean_bulk = bulk_corr(Rc, jc, nmax_c) * sgn_c
clean_edge = bulk_corr(Rc, 1, nmax_c) * sgn_c
rng = np.random.default_rng(3); dis = []
for r in range(40):
    Rr = R_matrix(J0 * rng.uniform(0.5, 1.5, Lc - 1), (np.pi / 2) * (1 - eps_c) + rng.uniform(-0.15, 0.15, Lc))
    dis.append(bulk_corr(Rr, jc, nmax_c) * sgn_c)
dis = np.mean(dis, axis=0)
print("clean bulk", clean_bulk[[10, 50, 200]], "clean edge", clean_edge[100:].mean(), "dis bulk", dis[100:].mean())

# ---------- Plot ----------
plt.rcParams.update({"font.size": 9, "font.family": "serif", "mathtext.fontset": "cm"})
fig, ax = plt.subplots(1, 3, figsize=(7.0, 2.35))
for k in range(qe.shape[1]):
    ax[0].plot(eps_list, qe[:, k] / np.pi, '.', ms=0.6, color='0.55')
# highlight modes pinned at pi
pin = np.abs(qe - np.pi) < 2e-3
ax[0].plot(np.repeat(eps_list[:, None], qe.shape[1], 1)[pin], qe[pin] / np.pi, '.', ms=2.5, color='C3')
ax[0].set_xlabel(r'$\epsilon$'); ax[0].set_ylabel(r'$|\theta_k|/\pi$')
ax[0].set_title(r'(a) Floquet spectrum, $L=40$', fontsize=9)
ax[0].set_ylim(-0.02, 1.05)

cols = ['C0', 'C1', 'C2', 'C3']
for (eps, (C, b)), c in zip(edge_data.items(), cols):
    sg = (-1) ** np.arange(len(C))
    ax[1].plot(C * sg, color=c, lw=0.9, label=rf'$\epsilon={eps}$')
    ax[1].axhline(b, color=c, ls='--', lw=0.7)
ax[1].set_xlabel(r'$n$'); ax[1].set_ylabel(r'$(-1)^n\,C_{\rm edge}(n)$')
ax[1].set_title(r'(b) edge spin, $\infty$-temperature', fontsize=9)
ax[1].legend(fontsize=7, frameon=False, loc='center right', bbox_to_anchor=(1.0,0.6)); ax[1].set_ylim(0, 1.05)

ax[2].plot(clean_edge, color='C3', lw=0.9, label='clean, edge')
ax[2].plot(clean_bulk, color='0.4', lw=0.9, label='clean, bulk')
ax[2].plot(dis, color='C0', lw=0.9, label='disordered, bulk')
ax[2].set_xlabel(r'$n$'); ax[2].set_ylabel(r'$(-1)^n\,C_j(n)$')
ax[2].set_title(r'(c) edge vs bulk, $\epsilon=0.1$', fontsize=9)
ax[2].legend(fontsize=7, frameon=False, loc='upper right'); ax[2].set_ylim(-0.05, 1.05)
plt.tight_layout(pad=0.3)
plt.savefig('fig_kicked_ising.pdf'); plt.savefig('fig_kicked_ising.png', dpi=150)
