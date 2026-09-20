"""O(j) exact stability decision for the path variety (tree_family.path(j)).

Group.  A = <reflection of hexagon i fixing r_i;  swap of f^c_1, f^c_2 in each P^2>, abelian, no
path flip.  Prop. 2.5 of notes/theory.md holds for any subgroup of Aut(Sigma) (the maximal
destabilizing subsheaf and the socle are invariant under every automorphism), so T_X is stable iff
the ray matroid is connected and rank(span S) - delta(S) > 0 for every A-invariant set S of rays
with span S != 0, N.  (Minimizing over all invariant S is the same as over invariant flats:
passing to the closure keeps the rank and adds weight.)

Structure.  N_Q = trivial isotypic part  <e_0, r_0, g_1, r_1, ..., g_{j-1}, r_{j-1}, e_j>
(g_c = f^c_1 + f^c_2)  plus the 2j-1 sign characters s_i = (1,2) in hexagon i and
h_c = f^c_1 - f^c_2.  The A-orbits of rays are {r_i}, {-r_i}, two pairs of other hexagon rays (each
pair spans the hexagon plane), {f^c_1, f^c_2}, {f^c_0}, and the four end rays.  So an invariant S
is, up to closure, a choice per site:
    hexagon i : nothing | the line (+-r_i) | the plane (six rays);
    P^2 no. c : {f_1,f_2} in or out (C)  x  f_0 in or out (Phi),  f^c_0 = -g_c + r_{c-1} - r_c;
    end 0     : e_0 in or out  x  f^0_0 = -e_0 - r_0 in or out;   end j likewise, f^j_0 = -e_j + r_{j-1}.
Rank = #planes + #C  (the sign characters)  +  rank of the chosen vectors in the trivial part
     = #chosen basis vectors (e_0, e_j, r_i, g_c) + #(f_0's whose own g_c or e is NOT chosen: a
       private coordinate) + rank of the remaining f_0's modulo the chosen basis vectors.
Those remaining f_0's are +-(r_{c-1} - r_c) with chosen r's set to 0: a graphic matroid on the
path r_0 - r_1 - ... with all chosen r's identified to one ground vertex.  Scanning left to right
with the state of r_c in {G (chosen), g (unchosen, joined to ground from the left), u (unchosen,
floating)}: an edge into a new unchosen vertex gains 1; an edge into ground gains 1 from u and 0
from G or g (it closes a cycle).
span S = N iff every hexagon is a plane, every P^2 has C, and each end has at least one ray.
Usage: path_flats.py [jmax]"""
import sys, os, time, random
from fractions import Fraction as Fr
from itertools import product
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, '..', '..', 'src'))
import path_transfer as PT

HEXOPT = (0, 1, 2)                                   # nothing, line, plane
BOOL2 = ((False, False), (True, False), (False, True), (True, True))

def weights(j):
    d, _ = PT.normalized_degrees(j)
    hexw = [(Fr(0), d[('hex', i, 'r')] + d[('hex', i, '-r')],
             d[('hex', i, 'r')] + d[('hex', i, '-r')] + 2*d[('hex', i, 'neg')] + 2*d[('hex', i, 'pos')]) for i in range(j)]
    basew = [d[('base', c)] for c in range(j+1)]      # weight of ONE ray of base c
    return hexw, basew

# ---------------------------------------------------------------- rank of one configuration (local rule)
def local_rank(j, end0, hexs, bases, endj):
    """end0 = (e_0 in, f^0_0 in); hexs[i] in {0,1,2}; bases[c-1] = (C, Phi) for c = 1..j-1; endj = (e_j in, f^j_0 in)."""
    rk = 0; e0, p0 = end0
    rk += e0; link = p0 and e0; rk += (p0 and not e0)
    state = None
    for c in range(j):
        if c >= 1:
            C, Phi = bases[c-1]; rk += 2*C; link = Phi and C; rk += (Phi and not C)
        if hexs[c] >= 1:
            rk += hexs[c]                              # r_c, and s_c for the plane
            if link and state == 'u': rk += 1
            state = 'G'
        else:
            if link: rk += 1; state = 'g' if (c == 0 or state in ('G', 'g')) else 'u'     # at c = 0 the edge comes from ground
            else: state = 'u'
    ej, pj = endj; rk += ej; rk += (pj and not ej)
    if pj and ej and state == 'u': rk += 1
    return rk

def config_rays(G, end0, hexs, bases, endj):
    idx = {l: q for q, l in enumerate(G.label)}; j = G.j; out = []
    if end0[0]: out.append(idx[('base', 0, 1)])
    if end0[1]: out.append(idx[('base', 0, 0)])
    for i, h in enumerate(hexs):
        if h == 1: out += [idx[('hex', i, (1, 0))], idx[('hex', i, (-1, 0))]]
        if h == 2: out += [idx[('hex', i, hh)] for hh in PT.KIND]
    for c, (C, Phi) in enumerate(bases, start=1):
        if C: out += [idx[('base', c, 1)], idx[('base', c, 2)]]
        if Phi: out.append(idx[('base', c, 0)])
    if endj[0]: out.append(idx[('base', j, 1)])
    if endj[1]: out.append(idx[('base', j, 0)])
    return out

# ---------------------------------------------------------------- the DP
def min_margin(j, hexw=None, basew=None):
    """Minimum of rank - delta over nonempty invariant S with span S != N, and a minimizing configuration.
    hexw[i] = (0, weight of the line, weight of the plane), basew[c] = weight of ONE ray of base c;
    default: the exact normalized degrees of the path with j hexagons."""
    if hexw is None: hexw, basew = weights(j)
    # DP state: (state of r_c in 'G','g','u', nonempty flag, full-type flag) -> (value, backpointer)
    cur = {}
    for end0 in BOOL2:
        e0, p0 = end0
        val = Fr(int(e0) + int(p0 and not e0)) - basew[0]*(int(e0) + int(p0))
        key = ('link' if (p0 and e0) else 'nolink', e0 or p0, e0 or p0)
        if key not in cur or val < cur[key][0]: cur[key] = (val, (end0,))
    def step_hex(cur, c):
        new = {}
        for (st, nonempty, full), (val, hist) in cur.items():
            link = st[0] == 'link' if isinstance(st, tuple) else st == 'link'
            prev = st[1] if isinstance(st, tuple) else None
            for h in HEXOPT:
                v = val - hexw[c][h]
                if h >= 1:
                    v += h + (1 if (link and prev == 'u') else 0); ns = 'G'
                else:
                    if link: v += 1; ns = 'g' if prev in (None, 'G', 'g') else 'u'
                    else: ns = 'u'
                key = (ns, nonempty or h >= 1, full and h == 2)
                if key not in new or v < new[key][0]: new[key] = (v, hist + (h,))
        return new
    def step_base(cur, c):
        new = {}
        for (st, nonempty, full), (val, hist) in cur.items():
            for C, Phi in BOOL2:
                v = val + 2*int(C) + int(Phi and not C) - basew[c]*(2*int(C) + int(Phi))
                key = (('link' if (Phi and C) else 'nolink', st), nonempty or C or Phi, full and C)
                if key not in new or v < new[key][0]: new[key] = (v, hist + ((C, Phi),))
        return new
    cur = step_hex(cur, 0)
    for c in range(1, j):
        cur = step_base(cur, c); cur = step_hex(cur, c)
    best = None
    for (st, nonempty, full), (val, hist) in cur.items():
        for ej, pj in BOOL2:
            v = val + int(ej) + int(pj and not ej) + (1 if (pj and ej and st == 'u') else 0) - basew[j]*(int(ej) + int(pj))
            if not (nonempty or ej or pj): continue                    # S empty
            if full and (ej or pj): continue                            # span S = N
            if best is None or v < best[0]: best = (v, hist + ((ej, pj),))
    return best

def describe(j, hist):
    end0, rest, endj = hist[0], hist[1:-1], hist[-1]
    hexs = [rest[0]] + [rest[2*c] for c in range(1, j)]; bases = [rest[2*c-1] for c in range(1, j)]
    parts = []
    if any(end0): parts.append("end0:" + ('e' if end0[0] else '') + ('f0' if end0[1] else ''))
    for i, h in enumerate(hexs):
        if i >= 1 and any(bases[i-1]): parts.append(f"P2_{i}:" + ('f12' if bases[i-1][0] else '') + ('f0' if bases[i-1][1] else ''))
        if h: parts.append(f"hex{i}:" + ('line' if h == 1 else 'plane'))
    if any(endj): parts.append(f"end{j}:" + ('e' if endj[0] else '') + ('f0' if endj[1] else ''))
    return ' '.join(parts), end0, hexs, bases, endj

# ---------------------------------------------------------------- Gate 3
def gate3(seed=1, n_random=20000):
    import tree_family as TF
    from toric_stability import rank as exact_rank
    rng = random.Random(seed)
    print("== Gate 3(i): local rank rule against exact rank of the chosen rays", flush=True)
    for j in range(1, 7):
        G = TF.path(j); t0 = time.time()
        exhaustive = 16 * 3**j * 4**(j-1) <= 7000
        if exhaustive: sample = list(product(BOOL2, product(HEXOPT, repeat=j), product(BOOL2, repeat=j-1), BOOL2))
        else: sample = [(rng.choice(BOOL2), tuple(rng.choice(HEXOPT) for _ in range(j)), tuple(rng.choice(BOOL2) for _ in range(j-1)), rng.choice(BOOL2))
                        for _ in range(n_random // 3)]
        bad = 0
        for end0, hexs, bases, endj in sample:
            idx = config_rays(G, end0, hexs, bases, endj)
            ex = exact_rank([list(G.rays[q]) for q in idx]) if idx else 0
            if ex != local_rank(j, end0, hexs, bases, endj): bad += 1
        print(f"  j = {j}: {len(sample)} configurations ({'all' if exhaustive else 'random'}), disagreements: {bad}  ({time.time()-t0:.0f}s)", flush=True)
        assert bad == 0
    print("== Gate 3(ii): DP minimum against the exhaustive minimum with exact ranks (all configurations)", flush=True)
    for j in range(1, 4):
        G = TF.path(j); hexw, basew = weights(j); best = None; n = 4*j
        for end0, hexs, bases, endj in product(BOOL2, product(HEXOPT, repeat=j), product(BOOL2, repeat=j-1), BOOL2):
            idx = config_rays(G, end0, hexs, bases, endj)
            if not idx: continue
            rk = exact_rank([list(G.rays[q]) for q in idx])
            if rk == n: continue
            w = basew[0]*sum(end0) + basew[j]*sum(endj) + sum(hexw[i][h] for i, h in enumerate(hexs)) + sum(basew[c]*(2*C + Phi) for c, (C, Phi) in enumerate(bases, start=1))
            if best is None or rk - w < best: best = rk - w
        dp = min_margin(j)[0]
        print(f"  j = {j}: exhaustive minimum {float(best):.8f}, DP minimum {float(dp):.8f}, equal: {best == dp}", flush=True)
        assert best == dp
    print("== Gate 3(iii): with the path flip the minimum is attained on both end caps: 2 x DP = tree_family minimum", flush=True)
    for j in range(2, 7):
        r = TF.path(j).verdict(); dp = min_margin(j)[0]
        print(f"  j = {j}: tree_family (flip on) {float(r['min_margin']):.8f};  2 x DP {float(2*dp):.8f};  equal: {r['min_margin'] == 2*dp}", flush=True)
        assert r['min_margin'] == 2*dp

if __name__ == '__main__':
    jmax = int(sys.argv[1]) if len(sys.argv) > 1 else 400
    gate3()
    print("\n== exact minimum margin over all A-invariant subspaces (no flip), by DP", flush=True)
    for j in list(range(1, 41)) + [50, 60, 80, 100, 150, 200, 300, 400]:
        if j > jmax: break
        t0 = time.time(); v, hist = min_margin(j); txt = describe(j, hist)[0]
        print(f"  j = {j:3d}  n = {4*j:4d}  rho = n+{j+1:<3d}  min margin = {float(v):.10f}  {'STABLE' if v > 0 else 'NOT STABLE'}   at [{txt}]  ({time.time()-t0:.1f}s)", flush=True)
