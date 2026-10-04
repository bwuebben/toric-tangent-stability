"""Exact stability certificates that use neither subspace enumeration nor the symmetry reduction.

Let E be the set of rays, r the rank function of the ray matroid (rank of the span), and
delta_v = n d_v / D the normalized degrees (delta > 0, delta(E) = n = r(E)).

Lemma A.  T_X is (-K)-stable  iff  delta(S) < r(S) for every proper nonempty S in E.
   (Klyachko / Hering-Nill-Suess: stable iff sum_{v in V} d_v < (dim V) D/n for every proper
   nonzero subspace V.  Given V, S = rays in V has r(S) <= dim V, and S = E would force V = N.
   Conversely apply the criterion to V = span S, which contains S, and V = N iff r(S) = n; if
   r(S) = n and S != E then delta(S) < delta(E) = n since delta > 0.)
Lemma B.  Suppose delta = sum_k lam_k 1_{B_k} with bases B_k and all lam_k > 0.  Then
   delta(S) <= r(S) for all S, with equality iff every B_k meets S in r(S) elements.  If moreover
   the digraph  e -> e'  (some k:  e not in B_k, e' in B_k, B_k - e' + e a basis)  is strongly
   connected, no proper nonempty S is tight: an arc e -> e' with e in S, e' not in S would give a
   basis meeting S in r(S) + 1 elements.  Hence T_X is stable.

Search is in floating point (column generation: master LP by HiGHS, pricing by the greedy
max-weight basis with ranks mod a 61-bit prime); the certificate is then re-derived and checked in
exact rational arithmetic, and src/verify_certificate.py re-checks the JSON from scratch.  When
delta is outside the base polytope the LP dual y separates it, and one of the level sets of y
violates delta(S) <= r(S); its closure is reported (exact ranks).
Usage: matroid_certificate.py [gate2 | name ...]"""
import sys, os, json, time, random
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import linprog
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
from toric_stability import rank as exact_rank, rref, in_span, mat_inverse
P = (1 << 61) - 1

# ---------------------------------------------------------------- greedy bases (ranks mod P)
def greedy_basis(rays, weights):
    """Indices of a maximum-weight basis (greedy; independence tested mod P)."""
    n = len(rays[0]); order = sorted(range(len(rays)), key=lambda q: -weights[q])
    rows, piv, basis = [], [], []
    for q in order:
        v = [x % P for x in rays[q]]
        for row, p in zip(rows, piv):
            if v[p]:
                f = v[p]; v = [(a - f*b) % P for a, b in zip(v, row)]
        p = next((k for k in range(n) if v[k]), None)
        if p is None: continue
        inv = pow(v[p], P-2, P); rows.append([(a*inv) % P for a in v]); piv.append(p); basis.append(q)
        if len(basis) == n: break
    return tuple(sorted(basis))

# ---------------------------------------------------------------- column generation
def decompose_float(rays, delta, seed=0, max_iter=4000, log=None):
    """Returns ('inside', pool, lam) or ('outside', pool, y)."""
    E = len(rays); n = len(rays[0]); rng = random.Random(seed); d = np.array([float(x) for x in delta])
    pool = []
    def add(B):
        if B not in pool: pool.append(B); return True
        return False
    add(greedy_basis(rays, list(d)))
    for _ in range(2*E): add(greedy_basis(rays, [rng.random() for _ in range(E)]))
    for it in range(max_iter):
        m = len(pool); A = np.zeros((E+1, m + 2*E))
        for k, B in enumerate(pool): A[list(B), k] = 1.0
        A[:E, m:m+E] = np.eye(E); A[:E, m+E:] = -np.eye(E); A[E, :m] = 1.0
        c = np.concatenate([np.zeros(m), np.ones(2*E)]); b = np.concatenate([d, [1.0]])
        res = linprog(c, A_eq=A, b_eq=b, bounds=(0, None), method='highs-ds')
        assert res.status == 0, res.message
        if res.fun < 1e-11: return 'inside', pool, res.x[:m]
        y = res.eqlin.marginals[:E]; sigma = res.eqlin.marginals[E]
        B = greedy_basis(rays, list(y))
        if sum(y[q] for q in B) + sigma <= 1e-9 or not add(B):
            return 'outside', pool, y
        if log and it % 50 == 0: log(f"    iteration {it}: residual {res.fun:.3e}, pool {len(pool)}")
    raise RuntimeError("column generation did not converge")

def exact_lambda(pool, support, delta):
    """Solve sum_k lam_k 1_{B_k} = delta exactly on the given support; None if inconsistent or not unique."""
    E = len(delta); s = len(support)
    M = [[Fr(1 if e in pool[k] else 0) for k in support] + [delta[e]] for e in range(E)]
    r = 0; piv = []
    for col in range(s):
        p = next((i for i in range(r, E) if M[i][col] != 0), None)
        if p is None: return None
        M[r], M[p] = M[p], M[r]; inv = 1 / M[r][col]; M[r] = [x*inv for x in M[r]]
        for i in range(E):
            if i != r and M[i][col] != 0:
                f = M[i][col]; M[i] = [a - f*b for a, b in zip(M[i], M[r])]
        piv.append(col); r += 1
    if any(M[i][s] != 0 for i in range(r, E)): return None
    return [M[i][s] for i in range(s)]

def exchange_digraph(rays, bases):
    """Arcs e -> e' : e not in B, e' in the fundamental circuit of e with respect to B (exact)."""
    E = len(rays); arcs = [set() for _ in range(E)]
    for B in bases:
        Binv = mat_inverse([list(rays[q]) for q in B])          # rows of B times Binv = identity
        inB = set(B)
        for e in range(E):
            if e in inB: continue
            coeff = [sum(Fr(rays[e][t]) * Binv[t][k] for t in range(len(B))) for k in range(len(B))]   # e = sum coeff_k B_k
            for k, cf in enumerate(coeff):
                if cf != 0: arcs[e].add(B[k])
    return arcs

def strongly_connected(arcs):
    E = len(arcs)
    def reach(adj):
        seen = {0}; stack = [0]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen: seen.add(y); stack.append(y)
        return len(seen) == E
    rev = [set() for _ in range(E)]
    for x in range(E):
        for y in arcs[x]: rev[y].add(x)
    return reach(arcs) and reach(rev)

def certify(rays, delta, name='', tries=12, log=print):
    """Returns a dict: status 'stable (certified)' with the certificate, or 'not certified' with a violated set."""
    E = len(rays); n = len(rays[0]); assert sum(delta) == n and all(x > 0 for x in delta)
    assert exact_rank([list(v) for v in rays]) == n
    acc = {}; n_dec = 0
    for t in range(tries):
        status, pool, out = decompose_float(rays, delta, seed=t, log=log if t == 0 else None)
        if status == 'outside':
            y = out; order = sorted(range(E), key=lambda q: -y[q]); worst = None
            for k in range(1, E):
                if y[order[k-1]] - y[order[k]] < 1e-12: continue
                S = order[:k]; rk, canon = rref([list(rays[q]) for q in S])
                if rk == n: continue
                closure = [q for q in range(E) if in_span(canon, rays[q])]
                viol = sum(delta[q] for q in closure) - rk
                if worst is None or viol > worst[0]: worst = (viol, rk, closure)
            return dict(name=name, status='not certified', violation=worst)
        support = [k for k, x in enumerate(out) if x > 1e-9]
        lam = exact_lambda(pool, support, delta)
        if lam is None or any(x <= 0 for x in lam): continue
        n_dec += 1
        for k, x in zip(support, lam): acc[pool[k]] = acc.get(pool[k], 0) + x
        bases = sorted(acc); lams = [acc[B] / n_dec for B in bases]
        for B in bases: assert exact_rank([list(rays[q]) for q in B]) == n
        assert all(sum(l for B, l in zip(bases, lams) if e in B) == delta[e] for e in range(E)) and sum(lams) == 1
        if strongly_connected(exchange_digraph(rays, bases)):
            return dict(name=name, status='stable (certified)', bases=bases, lambdas=lams, decompositions=n_dec)
    return dict(name=name, status='not certified', violation=None)

def save(cert, rays, delta, path, degree_source):
    json.dump({'name': cert['name'], 'degree_source': degree_source, 'rays': [list(v) for v in rays],
               'delta': [f"{x.numerator}/{x.denominator}" for x in delta],
               'bases': [list(B) for B in cert['bases']],
               'lambdas': [f"{x.numerator}/{x.denominator}" for x in cert['lambdas']]}, open(path, 'w'))

# ---------------------------------------------------------------- cases
def graph_case(G, use_transfer_j=None):
    if use_transfer_j:
        import path_transfer as PT
        d, _ = PT.normalized_degrees(use_transfer_j)
        delta = [d[('hex', l[1], PT.KIND[l[2]])] if l[0] == 'hex' else d[('base', l[1])] for l in G.label]
        return G.rays, delta, "path_transfer.normalized_degrees (3x3 transfer matrices; Gate 1: equal to tree_family for j <= 7)"
    deg, D = G.degrees()
    return G.rays, [G.n * x / D for x in deg], "tree_family.HexGraph.degrees (slice integrals; equal to the vertex formula in dims 8, 12, 16)"

def describe(G, idx):
    lab = [G.label[q] for q in idx]
    hx = sorted({l[1] for l in lab if l[0] == 'hex' and sum(1 for x in lab if x[:2] == l[:2]) == 6})
    bs = sorted({l[1] for l in lab if l[0] == 'base' and sum(1 for x in lab if x[:2] == l[:2]) == G.b[l[1]] + 1})
    other = [l for l in lab if not (l[0] == 'hex' and l[1] in hx) and not (l[0] == 'base' and l[1] in bs)]
    return f"hexagons {hx}, bases {bs}" + (f", other {other}" if other else '')

if __name__ == '__main__':
    import tree_family as T
    CASES = {'path2': (lambda: T.path(2), None), 'cycle2': (lambda: T.HexGraph([2, 2], [(0, 1), (1, 0)], "2-cycle on P^2, P^2"), None),
             'X14': (lambda: T.star(1, 4), None), 'X15': (lambda: T.star(1, 5), None), 'spider112': (lambda: T.spider(1, 1, 2), None),
             'path3': (lambda: T.path(3), None), 'path4': (lambda: T.path(4), None), 'X33': (lambda: T.star(3, 3), None),
             'spider222': (lambda: T.spider(2, 2, 2), None), 'path8': (lambda: T.path(8), 8),
             'spider122': (lambda: T.spider(1, 2, 2), None), 'spider113': (lambda: T.spider(1, 1, 3), None)}
    names = sys.argv[1:] or ['gate2']
    if names == ['gate2']: names = ['path2', 'cycle2', 'X14', 'X15', 'spider112']
    os.makedirs(os.path.join(HERE, '..', 'output', 'certificates'), exist_ok=True)
    for nm in names:
        if nm.startswith('path') and nm[4:].isdigit() and nm not in CASES:       # pathN: transfer-matrix degrees
            CASES[nm] = ((lambda N=int(nm[4:]): T.path(N)), int(nm[4:]))
        make, tj = CASES[nm]; G = make(); t0 = time.time()
        rays, delta, src = graph_case(G, tj)
        cert = certify(rays, delta, name=G.name, log=lambda s: None)
        if cert['status'].startswith('stable'):
            path = os.path.join(HERE, '..', 'output', 'certificates', nm + '.json'); save(cert, rays, delta, path, src)
            print(f"{G.name}: n = {G.n}, {len(rays)} rays -> STABLE, certified: {len(cert['bases'])} bases, min lambda = {float(min(cert['lambdas'])):.3e}, "
                  f"exchange digraph strongly connected; written to output/certificates/{nm}.json  ({time.time()-t0:.0f}s)", flush=True)
        else:
            v = cert['violation']
            msg = f"violated set: rank {v[1]}, delta - rank = {float(v[0]):+.5f}, = {describe(G, v[2])}" if v else "no decomposition with strongly connected exchange digraph found"
            print(f"{G.name}: n = {G.n}, {len(rays)} rays -> NOT certified; {msg}  ({time.time()-t0:.0f}s)", flush=True)
