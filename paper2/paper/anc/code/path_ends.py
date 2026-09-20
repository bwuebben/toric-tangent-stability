"""Paths whose two end vertices are P^{b0} and P^{bj}, b0, bj in {1,2,3} (inner vertices P^2):
X_j^{(b0,bj)} = HexGraph([b0, 2, ..., 2, bj], edges i -> i+1).  n = 4j - 2 + b0 + bj,
rho - n = j + 3 - b0 - bj.  The saturated path of theory.md section 11 is the type (1,1); the types
(2,1), (2,2) or (3,1), (3,2) give n = 1, 2, 3 mod 4.  Only the ends change:

  degrees   ell_0 = (b0+1-a)^{b0}/b0!, rho_0 = (bj+1+a)^{bj}/bj!  (degree <= 3, so 4x4 transfer
            matrices on 1, a, a^2, a^3; the kernel (3+a-a')^2/2 and the inner sites are unchanged);
            a ray of the end P^b has delta = int (b/L) ell_0 rho_{j-1} / int ell_0 rho_{j-1}, i.e.
            ell_0 replaced by (b0+1-a)^{b0-1}/(b0-1)!.
  group     A also contains S_b permuting f_1..f_b of an end P^b.  span(f_1..f_b) = <g> + standard
            representation (irreducible, multiplicity one in N), so an invariant set contains all of
            f_1..f_b or none: the end options are  split (rank b, weight b*delta)  x  twisted ray
            f_0 = -g -+ r (private coordinate if not split, else an edge ground - r).
  full      span S = N needs every hexagon a plane, every inner P^2 split, and each end split
            (for b = 1: split or twisted ray).
Everything else (contraction 2/7, bulk weights, period automaton, gamma_cyc, pigeonhole, budget) is
that of section 11.  This module re-implements the chain for general ends and checks itself
against tree_family (degrees), exact ranks (rank rule), exhaustive minima (DP), and the reviewed
modules path_flats / path_all_j at type (1,1).
Usage: path_ends.py            (all gates, then the all-j computation for every type)"""
import sys, os, time, random, math
from fractions import Fraction as Fr
from itertools import product
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
import path_transfer as PT
import path_flats as PF
import path_all_j as PA

N4 = 4
W = [[PT.I(p+q) for q in range(N4)] for p in range(N4)]
WH = {k: [[f(p+q) for q in range(N4)] for p in range(N4)] for k, f in PT.EDGE.items()}
pad = lambda M: [[M[p][q] if p < 3 and q < 3 else Fr(0) for q in range(N4)] for p in range(N4)]
KAP, KAP1 = pad(PT.KAP), pad(PT.KAP1)
mm = lambda A, B: [[sum(A[i][k]*B[k][j] for k in range(N4)) for j in range(N4)] for i in range(N4)]
vm = lambda v, A: [sum(v[k]*A[k][j] for k in range(N4)) for j in range(N4)]
mv = lambda A, v: [sum(A[i][k]*v[k] for k in range(N4)) for i in range(N4)]
dot = lambda u, v: sum(x*y for x, y in zip(u, v))
T, TR = mm(W, KAP), mm(KAP, W)

def end_vec(b, sign, drop=0):
    """Coefficients (in 1, a, a^2, a^3) of (b+1 + sign*a)^(b-drop) / (b-drop)!."""
    p = b - drop; out = [Fr(0)]*N4
    for q in range(p+1): out[q] = Fr(math.comb(p, q) * (b+1)**(p-q) * sign**q, math.factorial(p))
    return out

def chain_vectors(imax, b0, bj):
    L, R = [end_vec(b0, -1)], [end_vec(bj, +1)]
    for _ in range(imax): L.append(vm(L[-1], T)); R.append(mv(TR, R[-1]))
    return L, R

def hex_weights(Li, Rk):
    Z = dot(vm(Li, W), Rk); d = {k: dot(vm(Li, WH[k]), Rk) / Z for k in WH}
    line = d['r'] + d['-r']; return (Fr(0), line, line + 2*d['neg'] + 2*d['pos']), d

def base_weight(Li, Rk):
    return dot(vm(vm(vm(Li, W), KAP1), W), Rk) / dot(vm(vm(vm(Li, W), KAP), W), Rk)

def normalized_degrees(j, b0, bj):
    n = 4*j - 2 + b0 + bj; L, R = chain_vectors(j-1, b0, bj); d = {}
    Z = dot(vm(L[j-1], W), R[0])
    for i in range(j):
        _, dd = hex_weights(L[i], R[j-1-i])
        for k, x in dd.items(): d[('hex', i, k)] = x
    for c in range(1, j): d[('base', c)] = base_weight(L[c-1], R[j-1-c])
    d[('base', 0)] = dot(vm(end_vec(b0, -1, 1), W), R[j-1]) / dot(vm(L[0], W), R[j-1])
    d[('base', j)] = dot(vm(L[j-1], W), end_vec(bj, +1, 1)) / Z
    total = sum(d[('hex', i, 'r')] + d[('hex', i, '-r')] + 2*d[('hex', i, 'neg')] + 2*d[('hex', i, 'pos')] for i in range(j)) \
            + 3*sum(d[('base', c)] for c in range(1, j)) + (b0+1)*d[('base', 0)] + (bj+1)*d[('base', j)]
    assert total == n, (total, n)
    return d, math.factorial(n) * Z

def weights(j, b0, bj):
    d, _ = normalized_degrees(j, b0, bj)
    hexw = [(Fr(0), d[('hex', i, 'r')] + d[('hex', i, '-r')],
             d[('hex', i, 'r')] + d[('hex', i, '-r')] + 2*d[('hex', i, 'neg')] + 2*d[('hex', i, 'pos')]) for i in range(j)]
    return hexw, [d[('base', c)] for c in range(j+1)]

# ---------------------------------------------------------------- configurations
def local_rank(j, b0, bj, end0, hexs, bases, endj):
    """end = (split, twisted ray); hexs[i] in {0,1,2}; bases[c-1] = (split, twisted ray) for the inner P^2's."""
    sp0, tw0 = end0; rk = b0*sp0 + (tw0 and not sp0); link = tw0 and sp0; state = None
    for c in range(j):
        if c >= 1:
            C, Phi = bases[c-1]; rk += 2*C + (Phi and not C); link = Phi and C
        if hexs[c] >= 1:
            rk += hexs[c] + (1 if (link and state == 'u') else 0); state = 'G'
        elif link: rk += 1; state = 'g' if (c == 0 or state in ('G', 'g')) else 'u'
        else: state = 'u'
    spj, twj = endj; rk += bj*spj + (twj and not spj) + (1 if (twj and spj and state == 'u') else 0)
    return rk

def config_rays(G, b0, bj, end0, hexs, bases, endj):
    idx = {l: q for q, l in enumerate(G.label)}; j = G.j; out = []
    for c, b, (sp, tw) in ((0, b0, end0), (j, bj, endj)):
        if sp: out += [idx[('base', c, t)] for t in range(1, b+1)]
        if tw: out.append(idx[('base', c, 0)])
    for i, h in enumerate(hexs):
        if h == 1: out += [idx[('hex', i, (1, 0))], idx[('hex', i, (-1, 0))]]
        if h == 2: out += [idx[('hex', i, hh)] for hh in PT.KIND]
    for c, (C, Phi) in enumerate(bases, start=1):
        if C: out += [idx[('base', c, 1)], idx[('base', c, 2)]]
        if Phi: out.append(idx[('base', c, 0)])
    return out

def config_weight(j, b0, bj, hexw, basew, end0, hexs, bases, endj):
    return basew[0]*(b0*end0[0] + end0[1]) + basew[j]*(bj*endj[0] + endj[1]) + sum(hexw[i][h] for i, h in enumerate(hexs)) \
           + sum(basew[c]*(2*C + Phi) for c, (C, Phi) in enumerate(bases, start=1))

def min_margin(j, b0, bj, hexw=None, basew=None):
    if hexw is None: hexw, basew = weights(j, b0, bj)
    cur = {}
    for sp, tw in PF.BOOL2:
        val = Fr(b0*int(sp) + int(tw and not sp)) - basew[0]*(b0*int(sp) + int(tw))
        key = ('link' if (tw and sp) else 'nolink', sp or tw, sp or (b0 == 1 and tw))
        if key not in cur or val < cur[key][0]: cur[key] = (val, ((sp, tw),))
    def step_hex(cur, c):
        new = {}
        for (st, nonempty, full), (val, hist) in cur.items():
            link = (st[0] if isinstance(st, tuple) else st) == 'link'; prev = st[1] if isinstance(st, tuple) else None
            for h in PF.HEXOPT:
                v = val - hexw[c][h]
                if h >= 1: v += h + (1 if (link and prev == 'u') else 0); ns = 'G'
                elif link: v += 1; ns = 'g' if prev in (None, 'G', 'g') else 'u'
                else: ns = 'u'
                key = (ns, nonempty or h >= 1, full and h == 2)
                if key not in new or v < new[key][0]: new[key] = (v, hist + (h,))
        return new
    def step_base(cur, c):
        new = {}
        for (st, nonempty, full), (val, hist) in cur.items():
            for C, Phi in PF.BOOL2:
                v = val + 2*int(C) + int(Phi and not C) - basew[c]*(2*int(C) + int(Phi))
                key = (('link' if (Phi and C) else 'nolink', st), nonempty or C or Phi, full and C)
                if key not in new or v < new[key][0]: new[key] = (v, hist + ((C, Phi),))
        return new
    cur = step_hex(cur, 0)
    for c in range(1, j): cur = step_base(cur, c); cur = step_hex(cur, c)
    best = None
    for (st, nonempty, full), (val, hist) in cur.items():
        for sp, tw in PF.BOOL2:
            v = val + bj*int(sp) + int(tw and not sp) + (1 if (tw and sp and st == 'u') else 0) - basew[j]*(bj*int(sp) + int(tw))
            if not (nonempty or sp or tw): continue
            if full and (sp or (bj == 1 and tw)): continue
            if best is None or v < best[0]: best = (v, hist + ((sp, tw),))
    return best

# ---------------------------------------------------------------- gates
TYPES = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3)]

def gates():
    import tree_family as TF
    from toric_stability import rank as exact_rank
    rng = random.Random(7)
    print("== Gate E1: transfer-matrix degrees = tree_family degrees (exact), every ray, j <= 5", flush=True)
    for b0, bj in TYPES + [(1, 2), (1, 3), (2, 3)]:
        for j in range(1, 6):
            G = TF.HexGraph([b0] + [2]*(j-1) + [bj], [(i+1, i) for i in range(j)]); assert G.fano and G.n == 4*j - 2 + b0 + bj
            deg, D = G.degrees(); d, D2 = normalized_degrees(j, b0, bj)
            assert D == D2 and all(G.n*deg[q]/D == (d[('hex', l[1], PT.KIND[l[2]])] if l[0] == 'hex' else d[('base', l[1])]) for q, l in enumerate(G.label)), (b0, bj, j)
        print(f"   ends ({b0},{bj}): j = 1..5 agree", flush=True)
    print("== Gate E2/E3: rank rule = exact rank, and DP minimum = exhaustive minimum", flush=True)
    for b0, bj in TYPES:
        for j in range(1, 6):
            G = TF.HexGraph([b0] + [2]*(j-1) + [bj], [(i+1, i) for i in range(j)]); hexw, basew = weights(j, b0, bj)
            exhaustive = j <= 3
            space = product(PF.BOOL2, product(PF.HEXOPT, repeat=j), product(PF.BOOL2, repeat=j-1), PF.BOOL2) if exhaustive else \
                    ((rng.choice(PF.BOOL2), tuple(rng.choice(PF.HEXOPT) for _ in range(j)), tuple(rng.choice(PF.BOOL2) for _ in range(j-1)), rng.choice(PF.BOOL2)) for _ in range(2500))
            bad = 0; best = None; cnt = 0
            for end0, hexs, bases, endj in space:
                cnt += 1; idx = config_rays(G, b0, bj, end0, hexs, bases, endj)
                ex = exact_rank([list(G.rays[q]) for q in idx]) if idx else 0
                if ex != local_rank(j, b0, bj, end0, hexs, bases, endj): bad += 1
                if exhaustive and idx and ex < G.n:
                    m = ex - config_weight(j, b0, bj, hexw, basew, end0, hexs, bases, endj)
                    if best is None or m < best: best = m
            assert bad == 0, (b0, bj, j, bad)
            if exhaustive: assert best == min_margin(j, b0, bj)[0], (b0, bj, j)
        print(f"   ends ({b0},{bj}): rank rule exact on all configurations j <= 3 and 2,500 random ones for j = 4, 5; DP = exhaustive minimum j <= 3", flush=True)
    print("== Gate E4: type (1,1) reproduces the reviewed modules exactly", flush=True)
    for j in list(range(1, 31)) + [60, 84]:
        assert min_margin(j, 1, 1)[0] == PF.min_margin(j)[0]
        d1, D1 = normalized_degrees(j, 1, 1); d2, D2 = PT.normalized_degrees(j); assert D1 == D2 and d1 == d2
    print("   degrees and DP minima identical to path_transfer / path_flats for j = 1..30, 60, 84", flush=True)

# ---------------------------------------------------------------- all j
def all_j(b0, bj, q=40, KP=80):
    tau = Fr(2, 7); Delta = Fr(11757, 10000)
    L, R = chain_vectors(KP, b0, bj)
    hw, _ = hex_weights(L[KP], R[KP]); bw = base_weight(L[KP], R[KP])
    (gc, states, opts), _ = PA.gamma_cyc(PA.period_automaton(hw, bw)); assert gc > 0
    def model(m):
        j = 2*q + m; zone = lambda pos: 'left' if pos <= q-1 else ('deep' if pos <= q+m-1 else 'right')
        hexw = [hex_weights(L[i] if zone(i) == 'left' else L[KP], R[j-1-i] if zone(i) == 'right' else R[KP])[0] for i in range(j)]
        basew = [dot(vm(end_vec(b0, -1, 1), W), R[KP]) / dot(vm(L[0], W), R[KP])]
        basew += [base_weight(L[c-1] if zone(c) == 'left' else L[KP], R[j-1-c] if zone(c) == 'right' else R[KP]) for c in range(1, j)]
        basew.append(dot(vm(L[KP], W), end_vec(bj, +1, 1)) / dot(vm(L[KP], W), R[0]))
        assert all(h[1] < 2 and h[2] - h[1] < 4 for h in hexw) and all(0 < x < 1 for x in basew)
        return hexw, basew
    model_min, max_full, agree = [], Fr(0), Fr(0)
    for m in range(4):
        hexw, basew = model(m); j = 2*q + m; n = 4*j - 2 + b0 + bj
        v, hist = min_margin(j, b0, bj, hexw, basew); model_min.append(v)
        full = n - (sum(h[2] for h in hexw) + 3*sum(basew[1:j]) + (b0+1)*basew[0] + (bj+1)*basew[j]); max_full = max(max_full, abs(full))
        agree = max(agree, abs(v - min_margin(j, b0, bj)[0]))
    deep, zones, proxy = PA.uniform_error_terms(q, KP)
    budget = deep + zones + 2*proxy + max_full
    worst = None
    for j in range(1, 2*q + 5):
        v, hist = min_margin(j, b0, bj)
        dj, _ = normalized_degrees(j, b0, bj); assert all(0 < x < 1 for x in dj.values())
        if worst is None or v < worst[0]: worst = (v, j, hist)
    gstar = min([gc] + model_min)
    return dict(gamma_cyc=gc, model_min=min(model_min), agree=agree, budget=budget, worst=worst, bound=gstar - budget, max_full=max_full)

if __name__ == '__main__':
    gates()
    print("\n== all j: exact DP for every j <= 84, model paths M(0..3), error budget (q = 40, proxy index 80)", flush=True)
    print(" ends   n mod 4  rho-n   min over j<=84 (at j)        model minimum     |model - true| j=80..83   budget     bound for j>=80", flush=True)
    ok = True
    for b0, bj in TYPES:
        t0 = time.time(); r = all_j(b0, bj); s = (b0 + bj - 2) % 4
        stable_small = r['worst'][0] > 0; ok = ok and stable_small and r['bound'] > 0
        print(f" ({b0},{bj})    {s}      j{3-b0-bj:+d}    {float(r['worst'][0]):.12f} (j={r['worst'][1]:2d})   {float(r['model_min']):.12f}   {float(r['agree']):.1e}"
              f"                 {float(r['budget']):.1e}   {float(r['bound']):.12f}   {'STABLE for every j >= 1' if stable_small and r['bound'] > 0 else 'FAILS'}  ({time.time()-t0:.0f}s)", flush=True)
    print("\nCONCLUSION:", "every type is stable for every j >= 1." if ok else "some type fails — see the table.", flush=True)
