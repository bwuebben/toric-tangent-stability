"""The finite computations behind: the path variety is stable for every j  (notes/theory.md, section 11).

Notation.  ell_i (function of a_i) and rho_k (function of a_{j-1-k}) are the left and right chain
functions, ell_0 = 2 - a, rho_0 = 2 + a, each obtained from the previous one by the positive
integral operator with kernel K(a,a') = (3 + a - a')^2/2 against (2-|a|)da; as coefficient
vectors, L_i = lam^T (W kappa)^i and R_k = (kappa W)^k rho (path_transfer.py).  A normalized degree
at a site with parameters (i, k) is a ratio of two integrals of ell_i * rho_k against positive
measures:   hexagon i:  (i, j-1-i);   P^2 no. c:  (c-1, j-1-c);   end P^1's:  k = j-1, resp. i = j-1.

(H) Hilbert-metric contraction.  The cross-ratio K(a,a')K(b,b') / (K(a,b')K(b,a')) is at most
    81/25 on [-1,1]^4, so (Birkhoff-Hopf) the operator has projective diameter <= Delta = log(81/25)
    and contraction ratio tau = tanh(Delta/4) = 2/7.  Hence for i, i' >= 1 (incl. i' = infinity, the
    Perron eigenfunction) d_H(ell_i, ell_i') <= Delta tau^(min(i,i')-1), and a normalized degree
    changes by a factor within exp(+-(d_H of the ell's + d_H of the rho's)).
(Z) With bulk weights delta(inf,inf) the all-full period costs exactly 0 (sum delta = 4j for all j
    plus (H)); the all-empty period costs 0.
(C) Every other simple cycle of the period automaton on {G, g, u} costs >= gamma_cyc > 0.
(M) Model path M(m): K-site zones at both ends with weights delta(i,inf), delta(inf,k), and m bulk
    periods.  Removing closed walks between equal states (pigeonhole, m >= 4) reduces to m <= 3.
(E) Error budget between the true weights of the path with j = 2K + m hexagons and M(m), and
    between the reference weights (algebraic) and their rational proxies with index KP.
Conclusion: min margin(j) >= min(gamma_cyc, min_{m<=3} min margin M(m)) - budget for j >= 2K; and
the exact DP covers every j <= 2K + 4.
Usage: path_all_j.py [K] [KP]"""
import sys, os, time
from fractions import Fraction as Fr
from itertools import product, permutations
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import path_transfer as PT
import path_flats as PF

def chain_vectors(imax):
    L, R = [list(PT.LAM)], [list(PT.RHO)]
    for _ in range(imax): L.append(PT.vm(L[-1], PT.T)); R.append(PT.mv(PT.TR, R[-1]))
    return L, R

def hex_weights(Li, Rk):
    Z = PT.dot(PT.vm(Li, PT.W), Rk); d = {k: PT.dot(PT.vm(Li, PT.WH[k]), Rk) / Z for k in PT.WH}
    line = d['r'] + d['-r']; return (Fr(0), line, line + 2*d['neg'] + 2*d['pos'])

def base_weight(Li, Rk):            # one ray of a P^2 between ell_i and rho_k
    num = PT.dot(PT.vm(PT.vm(PT.vm(Li, PT.W), PT.KAP1), PT.W), Rk); den = PT.dot(PT.vm(PT.vm(PT.vm(Li, PT.W), PT.KAP), PT.W), Rk)
    return num / den

def model_weights(m, K, L, R, KP):
    """Weights of M(m), j = 2K + m hexagons, exactly as in the proof of Theorem 11.5; index 'infinity' -> proxy KP.
    Left zone: end cap and periods 1..K-1 (hexagons i <= K-1, P^2's c <= K-1): weights delta(i, inf), delta(c-1, inf).
    Deep periods p = K..K+m-1: BOTH the P^2 no. p and the hexagon p carry the bulk weights delta(inf, inf).
    Right zone: periods K+m..j-1 and the last end: delta(inf, j-1-i), delta(inf, j-1-c).
    (Zone membership is by position, not by the site's parameters: the P^2 no. K+m has left parameter K+m-1.)"""
    j = 2*K + m; hexw, basew = [], []
    def params(pos):                                   # pos = index of the period (hexagon i, or P^2 no. c)
        if pos <= K-1: return 'left'
        return 'deep' if pos <= K+m-1 else 'right'
    for i in range(j):
        z = params(i); hexw.append(hex_weights(L[i] if z == 'left' else L[KP], R[j-1-i] if z == 'right' else R[KP]))
    basew.append(PT.dot(PT.vm(PT.ONE, PT.W), R[KP]) / PT.dot(PT.vm(PT.LAM, PT.W), R[KP]))          # end 0: (2-a) -> 1, k = infinity
    for c in range(1, j):
        z = params(c); basew.append(base_weight(L[c-1] if z == 'left' else L[KP], R[j-1-c] if z == 'right' else R[KP]))
    basew.append(PT.dot(PT.vm(L[KP], PT.W), PT.ONE) / PT.dot(PT.vm(L[KP], PT.W), PT.RHO))          # end j: i = infinity
    assert all(h[1] < 2 and h[2] - h[1] < 4 for h in hexw) and all(0 < b < 1 for b in basew)       # aggregate finite-proxy diagnostics only
    return hexw, basew

def period_automaton(hw, bw):
    """Transitions of one bulk period [P^2, hexagon] between states of r: cost and next state."""
    trans = []
    for st in 'Ggu':
        for (C, Phi), h in product(PF.BOOL2, PF.HEXOPT):
            link = Phi and C; cost = 2*int(C) + int(Phi and not C) - bw*(2*int(C) + int(Phi)) - hw[h]
            if h >= 1: cost += h + (1 if (link and st == 'u') else 0); ns = 'G'
            elif link: cost += 1; ns = 'g' if st in 'Gg' else 'u'
            else: ns = 'u'
            trans.append((st, ns, (C, Phi, h), cost))
    return trans

def gamma_cyc(trans):
    zero = {('G', (True, True, 2)), ('u', (False, False, 0))}
    best = {}
    for st, ns, opt, cost in trans:
        if (st, opt) in zero: continue
        if (st, ns) not in best or cost < best[(st, ns)][0]: best[(st, ns)] = (cost, opt)
    cyc = []
    for r in (1, 2, 3):
        for states in permutations('Ggu', r):
            edges = [(states[t], states[(t+1) % r]) for t in range(r)]
            if all(e in best for e in edges): cyc.append((sum(best[e][0] for e in edges), states, [best[e][1] for e in edges]))
    return min(cyc, key=lambda c: c[0]), {k: v for k, v in best.items()}

def upper_exp_minus_1(x):           # e^x - 1 <= x + x^2 for 0 <= x <= 1
    assert 0 <= x <= 1; return x + x*x

def uniform_error_terms(q=40, proxy_index=80):
    """Rational bounds valid for all six endpoint types and all path lengths.

    Every chain function has max/min <= 9, including the initial end polynomials.
    A hexagon edge has mass 1 and the interior measure mass 3, so every normalized
    hexagon-ray weight is <= 9^2/3 = 27. Base-ray weights are <= 2. An end zone
    has at most 9*q+4 rays. No unproved all-length bound delta < 1 is used.
    """
    import math
    tau, Delta, B = Fr(2, 7), Fr(11757, 10000), Fr(27)
    # exp(Delta) > 81/25 proves Delta > log(81/25), using rationals only.
    assert sum(Delta**r / math.factorial(r) for r in range(8)) > Fr(81, 25)
    x, y = Delta*tau**(q-2), 2*Delta*tau**(proxy_index-1)
    ssum = 2*x / (1-tau)
    deep = B*9*ssum*(1+2*x)
    zones = B*2*(9*q+4)*upper_exp_minus_1(x)
    proxy = B*(9*(2*q+4)+4)*upper_exp_minus_1(y)
    return deep, zones, proxy

if __name__ == '__main__':
    K = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    KP = int(sys.argv[2]) if len(sys.argv) > 2 else 2*K
    tau = Fr(2, 7); Delta = Fr(11757, 10000)          # log(81/25) = 1.17557... < 1.1757
    uniform_error_terms(K, KP)  # includes a rational verification of Delta
    print(f"K = {K}, proxy index KP = {KP}, tau = {tau}, Delta <= {float(Delta)}", flush=True)

    print("\n== (H) cross-ratio of L = 3 + a - a' on a grid (the proof is by monotonicity; this is a check)", flush=True)
    g = [Fr(t, 20) - 1 for t in range(41)]; Lf = lambda a, b: 3 + a - b
    mx = max(Lf(a, a1)*Lf(b, b1) / (Lf(a, b1)*Lf(b, a1)) for a in g for b in g for a1 in (g[0], g[-1], g[20]) for b1 in (g[0], g[-1], g[20]))
    print(f"   max over the grid = {mx} (claimed supremum 9/5, squared 81/25)", flush=True); assert mx == Fr(9, 5)

    L, R = chain_vectors(KP)
    print("\n== (Z) bulk period weight at the proxy index (exactly 4 in the limit)", flush=True)
    hw = hex_weights(L[KP], R[KP]); bw = base_weight(L[KP], R[KP])
    print(f"   3 * (P^2 ray) + hexagon plane = 4 + {float(3*bw + hw[2] - 4):.3e}", flush=True)

    print("\n== (C) simple cycles of the period automaton, bulk weights", flush=True)
    (gc, states, opts), best = gamma_cyc(period_automaton(hw, bw))
    for k2 in sorted(best): print(f"   cheapest {k2[0]} -> {k2[1]} (zero loops excluded): {float(best[k2][0]):.6f}  option (C, Phi, hexagon) = {best[k2][1]}")
    print(f"   gamma_cyc = {float(gc):.8f}, on the cycle {states} with options {opts}", flush=True); assert gc > 0

    print("\n== (M) model paths M(m), m = 0..3: minimum over nontrivial invariant configurations", flush=True)
    model_min = []; max_full = Fr(0)
    for m in range(4):
        hexw, basew = model_weights(m, K, L, R, KP); j = 2*K + m
        v, hist = PF.min_margin(j, hexw, basew); model_min.append(v)
        full_cost = 4*j - (sum(h[2] for h in hexw) + 3*sum(basew[1:j]) + 2*basew[0] + 2*basew[j])
        tv, _ = PF.min_margin(j); max_full = max(max_full, abs(full_cost))
        print(f"   m = {m}: model minimum {float(v):.12f} at [{PF.describe(j, hist)[0]}];  all-full configuration costs {float(full_cost):+.3e};"
              f"  true path j = {j}: {float(tv):.12f}  (difference {float(v - tv):+.2e})", flush=True)
    gstar = min([gc] + model_min)

    print("\n== (E) error budget (exact rationals; uniform ray-weight bound 27)", flush=True)
    deep, zones, proxy = uniform_error_terms(K, KP)
    budget = deep + zones + 2*proxy + max_full
    print(f"   deep bulk (geometric tails) <= {float(deep):.3e};  end zones <= {float(zones):.3e};  proxies <= {float(proxy):.3e};"
          f"  |all-full cost on the model| <= {float(max_full):.3e}")
    print(f"   total budget <= {float(budget):.3e}", flush=True)

    print(f"\n== exact DP for every j <= 2K + 4 = {2*K+4}", flush=True); t0 = time.time(); worst = None
    for j in range(1, 2*K + 5):
        v, _ = PF.min_margin(j); assert v > 0
        dj, _ = PT.normalized_degrees(j); assert all(0 < x < 1 for x in dj.values())          # finite-path diagnostic; the uniform proof uses the bound 27
        if worst is None or v < worst[0]: worst = (v, j)
    print(f"   all stable; smallest minimum margin {float(worst[0]):.12f} at j = {worst[1]}  ({time.time()-t0:.0f}s)", flush=True)

    bound = gstar - budget; assert bound > 0
    print(f"\nCONCLUSION: for every j >= {2*K}: minimum margin >= min(gamma_cyc, model minima) - budget = {float(gstar):.12f} - {float(budget):.1e} > 0;"
          f" for j <= {2*K+4}: exact.  The path variety is stable for every j >= 1.", flush=True)
