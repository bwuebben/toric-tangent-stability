"""Exact degrees of the path variety for every j, by 3x3 transfer matrices.

The path with j hexagons (tree_family.path(j)): base factors c = 0..j (P^1 at the two ends, P^2
inside), hexagon i between base i and base i+1, a_i = <m_i, r_i> in [-1,1].  The slice of the
moment polytope over (m_i) has volume

    (2 - a_0) * prod_{c=1}^{j-1} K(a_{c-1}, a_c) * (2 + a_{j-1}),      K(a,a') = (3 + a - a')^2 / 2,

integrated against prod (2 - |a_i|) da_i.  K is a polynomial of degree 2 in each variable, so with
W = (int a^{p+q} (2-|a|) da)_{p,q<=2} and K = sum kappa[p][q] a^p a'^q,

    Z_j = lam^T (W kappa)^{j-1} W rho,    lam = (2,-1,0),  rho = (2,1,0),    D = (4j)! Z_j.

Normalized degrees delta_v = n d_v / D (they sum to n = 4j; T_X is stable iff delta(F) < rank F
for every proper nonempty flat F of the ray matroid) are ratios X_v / Z_j, X_v the same chain with
one factor replaced:
    ray h of hexagon i : the i-th W      -> W_h   (moments of the edge measure of h)
    a ray of P^2 no. c : the c-th kappa  -> kappa' (coefficients of 3 + a - a';  L^2/2 * 2/L = L)
    a ray of an end P^1: lam or rho      -> (1,0,0)
Usage: path_transfer.py [jmax]   (Gate 1 against tree_family for j <= 7, then the table to jmax)"""
import sys, os, math
from fractions import Fraction as Fr
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

def I(q): return Fr(0) if q % 2 else 2*(Fr(2, q+1) - Fr(1, q+2))          # int a^q (2-|a|) da
EDGE = {'r':   lambda q: Fr((-1)**q),        # a = -1
        '-r':  lambda q: Fr(1),              # a = +1
        'neg': lambda q: Fr((-1)**q, q+1),   # uniform on [-1,0]
        'pos': lambda q: Fr(1, q+1)}         # uniform on [0,1]
KIND = {(1,0): 'r', (-1,0): '-r', (1,1): 'neg', (0,-1): 'neg', (0,1): 'pos', (-1,-1): 'pos'}
W   = [[I(p+q) for q in range(3)] for p in range(3)]
WH  = {k: [[f(p+q) for q in range(3)] for p in range(3)] for k, f in EDGE.items()}
KAP  = [[Fr(9,2), Fr(-3), Fr(1,2)], [Fr(3), Fr(-1), Fr(0)], [Fr(1,2), Fr(0), Fr(0)]]   # (3+a-a')^2/2
KAP1 = [[Fr(3), Fr(-1), Fr(0)], [Fr(1), Fr(0), Fr(0)], [Fr(0), Fr(0), Fr(0)]]           # 3+a-a'
LAM, RHO, ONE = [Fr(2), Fr(-1), Fr(0)], [Fr(2), Fr(1), Fr(0)], [Fr(1), Fr(0), Fr(0)]

def mm(A, B): return [[sum(A[i][k]*B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
def vm(v, A): return [sum(v[k]*A[k][j] for k in range(3)) for j in range(3)]          # row vector * matrix
def mv(A, v): return [sum(A[i][k]*v[k] for k in range(3)) for i in range(3)]          # matrix * column
def dot(u, v): return sum(x*y for x, y in zip(u, v))
T, TR = mm(W, KAP), mm(KAP, W)                                                        # W kappa,  kappa W

def chains(j, lam=LAM, rho=RHO):
    """left[i] = lam^T (W kappa)^i,  right[k] = (kappa W)^k rho,  0 <= i,k <= j-1."""
    left, right = [list(lam)], [list(rho)]
    for _ in range(j-1): left.append(vm(left[-1], T)); right.append(mv(TR, right[-1]))
    return left, right

def normalized_degrees(j):
    """dict: ('hex', i, kind) -> delta of one ray of that kind;  ('base', c) -> delta of one ray of base c."""
    n = 4*j; left, right = chains(j)
    Z = dot(vm(left[j-1], W), RHO)
    d = {}
    for i in range(j):
        for k in WH: d[('hex', i, k)] = dot(vm(left[i], WH[k]), right[j-1-i]) / Z
    for c in range(1, j):
        d[('base', c)] = dot(vm(vm(vm(left[c-1], W), KAP1), W), right[j-1-c]) / Z
    l1, _ = chains(j, lam=ONE)
    d[('base', 0)] = dot(vm(l1[j-1], W), RHO) / Z            # (2 - a_0) -> 1
    d[('base', j)] = dot(vm(left[j-1], W), ONE) / Z          # (2 + a_{j-1}) -> 1
    total = sum(d[('hex', i, 'r')] + d[('hex', i, '-r')] + 2*d[('hex', i, 'neg')] + 2*d[('hex', i, 'pos')] for i in range(j)) \
            + 3*sum(d[('base', c)] for c in range(1, j)) + 2*d[('base', 0)] + 2*d[('base', j)]
    assert total == n, (total, n)
    return d, math.factorial(n) * Z

def gate1(jmax=7):
    import tree_family as TF
    for j in range(1, jmax+1):
        G = TF.path(j); deg, D = G.degrees(); d, D2 = normalized_degrees(j)
        ok = D == D2
        for q, l in enumerate(G.label):
            key = ('hex', l[1], KIND[l[2]]) if l[0] == 'hex' else ('base', l[1])
            ok = ok and (G.n * deg[q] / D == d[key])
        print(f"  j = {j}: D = {D2}; all {len(deg)} normalized degrees equal those of tree_family: {ok}", flush=True)
        assert ok

def bulk():
    """Perron data of T = W kappa (floating point, for orientation) and the exact characteristic polynomial."""
    import numpy as np
    tr = T[0][0] + T[1][1] + T[2][2]
    m2 = sum(T[i][i]*T[k][k] - T[i][k]*T[k][i] for i in range(3) for k in range(i+1, 3))
    dt = (T[0][0]*(T[1][1]*T[2][2]-T[1][2]*T[2][1]) - T[0][1]*(T[1][0]*T[2][2]-T[1][2]*T[2][0]) + T[0][2]*(T[1][0]*T[2][1]-T[1][1]*T[2][0]))
    print(f"  characteristic polynomial of T = W kappa:  x^3 - ({tr}) x^2 + ({m2}) x - ({dt})")
    A = np.array([[float(x) for x in row] for row in T]); ev, R = np.linalg.eig(A); evl, L = np.linalg.eig(A.T)
    order = np.argsort(-np.abs(ev)); print("  eigenvalues:", [complex(ev[o]) if abs(ev[o].imag) > 1e-12 else float(ev[o].real) for o in order],
                                          f"  |lambda_2/lambda_1| = {abs(ev[order[1]])/abs(ev[order[0]]):.6f}")
    # bulk limits: u^T T = Lam u^T (left), (kappa W) v = Lam v (right)
    u = np.real(L[:, np.argmax(np.abs(evl))]); B = np.array([[float(x) for x in row] for row in TR]); evr, Rr = np.linalg.eig(B)
    v = np.real(Rr[:, np.argmax(np.abs(evr))]); Lam = float(np.max(np.abs(ev)))
    Wf = np.array([[float(x) for x in row] for row in W]); K1 = np.array([[float(x) for x in row] for row in KAP1])
    den = u @ Wf @ v; lim = {k: float(u @ np.array([[float(x) for x in row] for row in WH[k]]) @ v / den) for k in WH}
    lim['P2'] = float(u @ Wf @ K1 @ Wf @ v / (Lam * den))
    per = lim['r'] + lim['-r'] + 2*lim['neg'] + 2*lim['pos'] + 3*lim['P2']
    print("  bulk limits of the normalized degrees:", {k: round(x, 8) for k, x in lim.items()}, f"  one period (6 hexagon rays + 3 P^2 rays) sums to {per:.8f} (dimension 4)")
    return lim

if __name__ == '__main__':
    jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    print("== Gate 1: transfer-matrix degrees against tree_family (polynomial expansion), exact", flush=True)
    gate1(7)
    print("\n== spectrum and bulk limits", flush=True); lim = bulk()
    print("\n== normalized degrees along the path (floats of exact rationals)", flush=True)
    print("   j | end P^1 | end hexagon: r, -r, neg, pos | first P^2 | middle hexagon: r, -r, neg, pos | middle P^2")
    for j in [2, 3, 4, 6, 8, 12, 20, 50, 100, 200, jmax]:
        if j > jmax: continue
        d, _ = normalized_degrees(j); mid = j // 2; f = lambda x: f"{float(x):.7f}"
        print(f"{j:4d} | {f(d[('base', 0)])} | " + ' '.join(f(d[('hex', 0, k)]) for k in ('r', '-r', 'neg', 'pos')) + f" | {f(d[('base', 1)])} | "
              + ' '.join(f(d[('hex', mid, k)]) for k in ('r', '-r', 'neg', 'pos')) + f" | {f(d[('base', max(1, mid))])}", flush=True)
