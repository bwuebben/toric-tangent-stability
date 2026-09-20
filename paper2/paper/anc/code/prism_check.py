#!/usr/bin/env python3
"""Exact validation of the PRISM STRUCTURE of the root-twist family.

Claim (Lemma, to enter the paper): the anticanonical polytope of
X(t_1..t_k) is the generalized prism

    P = { (x,y) : x in Q,  -1 - <x,t_i> <= y_i <= 1 },

Q = the dP3 hexagon {|x1|<=1, |x2|<=1, |x1+x2|<=1} in the M-lattice of the
fibre plane.  Writing l_t(x) = <x,t> and F(x) = prod_i (2 + l_{t_i}(x)):

    (-K)^n           = n!  * I_Q[F]
    deg D_{f_i+t_i}  = deg D_{-f_i} = (n-1)! * I_Q[F/(2+l_{t_i})]
    deg D_rho        = (n-1)! * I_{E_rho}[F]     (edge, lattice measure)
    mu(T_X)          = (n-1)! * I_Q[F]

and the divergence identity  I_{dQ}[F] = 2 I_Q[F] + I_Q[<x, grad F>]
turns every stability comparison into a hexagon integral.  This script
verifies ALL of these exactly (Fraction arithmetic) against the full-dimensional
harness (analyze_fast == analyze on the reference-sized cases) over:

  * the 8 canonical cases stored in output/root_twist.json (dims 4-8),
  * every k=1 and k=2 twist tuple over {roots} u {0} (dims 3-4),
  * random k=3 and k=4 tuples (dims 5-6),

and additionally checks, case by case:
  (a) all 6+2k facet degrees match the integral formulas exactly;
  (b) the REDUCED subspace list (W in {0, root line, N_f} + pair/single
      lifts) reproduces the harness max subsheaf slope exactly -- this is
      the "finitely many combinatorial types" classification;
  (c) mu(F_{N_f}) - mu(T) = ((n-1)!/2) * Phi with
      Phi = I_Q[ F * sum_i l_{t_i}/(2+l_{t_i}) ]   (divergence identity);
  (d) the full barycenter agrees with the prism formula, and hence
      KE  <=>  I_Q[F*x1] = I_Q[F*x2] = 0             (moment criterion).
"""

import json
import os
import random
import sys
from fractions import Fraction as Fr
from itertools import combinations, product as iproduct
from math import comb as binom, factorial

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from root_twist import root_twist_rays                      # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

R6 = [(1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1)]
KEY = {(1, 0): "a", (-1, 0): "-a", (0, 1): "b", (0, -1): "-b",
       (1, 1): "a+b", (-1, -1): "-(a+b)", (0, 0): "0"}
VEC = {v: k for k, v in KEY.items()}

# ---------------- exact polynomials in (x1,x2): dict[(a,b)] -> Fraction
def pmul(P, Q):
    R = {}
    for (a, b), c in P.items():
        for (d, e), f in Q.items():
            k = (a + d, b + e)
            R[k] = R.get(k, Fr(0)) + c * f
    return {k: v for k, v in R.items() if v}


def two_plus_ell(t):
    P = {(0, 0): Fr(2)}
    if t[0]:
        P[(1, 0)] = Fr(t[0])
    if t[1]:
        P[(0, 1)] = P.get((0, 1), Fr(0)) + Fr(t[1])
    return {k: v for k, v in P.items() if v}


def ell_poly(t):
    P = {}
    if t[0]:
        P[(1, 0)] = Fr(t[0])
    if t[1]:
        P[(0, 1)] = Fr(t[1])
    return P


# ---------------- exact monomial integrals over the hexagon
def _int_square(a, b):
    fa = Fr(2, a + 1) if a % 2 == 0 else Fr(0)
    fb = Fr(2, b + 1) if b % 2 == 0 else Fr(0)
    return fa * fb


def _int_Tplus(a, b):
    # triangle {0<=x<=1, 1-x<=y<=1} (the corner of the square cut by Q)
    beta = Fr(factorial(a) * factorial(b + 1), factorial(a + b + 2))
    return (Fr(1, a + 1) - beta) / (b + 1)


def _int_hex_mono(a, b):
    tp = _int_Tplus(a, b)
    tm = tp if (a + b) % 2 == 0 else -tp
    return _int_square(a, b) - tp - tm


def int_hex(P):
    return sum(c * _int_hex_mono(a, b) for (a, b), c in P.items())


# ---------------- exact edge integrals (lattice measure, param s in [0,1])
EDGES = {(1, 0): ((-1, 0), (0, 1)),     # E_a: x1=-1, x2 in [0,1]
         (-1, 0): ((1, -1), (0, 1)),    # E_-a: x1=1, x2 in [-1,0]
         (0, 1): ((0, -1), (1, 0)),     # E_b: x2=-1
         (0, -1): ((-1, 1), (1, 0)),    # E_-b: x2=1
         (1, 1): ((-1, 0), (1, -1)),    # E_{a+b}: x1+x2=-1
         (-1, -1): ((0, 1), (1, -1))}   # E_{-(a+b)}: x1+x2=1


def int_edge(P, rho):
    (p1, p2), (d1, d2) = EDGES[rho]
    coef = {}
    for (a, b), c in P.items():
        pa = [Fr(binom(a, j)) * Fr(p1) ** (a - j) * Fr(d1) ** j
              for j in range(a + 1)]
        pb = [Fr(binom(b, j)) * Fr(p2) ** (b - j) * Fr(d2) ** j
              for j in range(b + 1)]
        for i, ca in enumerate(pa):
            if not ca:
                continue
            for j, cb in enumerate(pb):
                if not cb:
                    continue
                coef[i + j] = coef.get(i + j, Fr(0)) + c * ca * cb
    return sum(c / (m + 1) for m, c in coef.items())


# ---------------- degree formulas from the prism lemma
def prism_degrees(twists):
    k = len(twists)
    n = k + 2
    F = {(0, 0): Fr(1)}
    for t in twists:
        F = pmul(F, two_plus_ell(t))
    intF = int_hex(F)
    total = factorial(n) * intF
    mu = Fr(factorial(n - 1)) * intF
    a = {rho: factorial(n - 1) * int_edge(F, rho) for rho in R6}
    b = []
    for i in range(k):
        Fi = {(0, 0): Fr(1)}
        for j, s in enumerate(twists):
            if j != i:
                Fi = pmul(Fi, two_plus_ell(s))
        b.append(factorial(n - 1) * int_hex(Fi))
    return a, b, total, mu, F


# ---------------- reduced subspace list (the classification)
def reduced_max_slope(twists, a, b):
    """Max slope over: W in {0, root line, N_f}, plus per-factor additions
    (skip / single lift / pair).  Pair for nonzero t_i needs t_i in W."""
    k = len(twists)
    n = k + 2
    lines = [((1, 0), (-1, 0)), ((0, 1), (0, -1)), ((1, 1), (-1, -1))]
    bases = [(Fr(0), 0, None)]
    for rp, rm in lines:
        bases.append((a[rp] + a[rm], 1, {rp, rm}))
    bases.append((sum(a.values()), 2, "all"))
    best = Fr(-1)
    for bnum, bdim, wroots in bases:
        opts = []
        for i, t in enumerate(twists):
            o = [(Fr(0), 0), (b[i], 1)]                    # skip / single
            pair_ok = (t == (0, 0)) or (wroots == "all") \
                or (wroots is not None and t in wroots)
            if pair_ok:
                o.append((2 * b[i], 1))                    # pair
            opts.append(o)
        for choice in iproduct(*opts):
            num = bnum + sum(c[0] for c in choice)
            dim = bdim + sum(c[1] for c in choice)
            if 1 <= dim <= n - 1:
                s = Fr(num, dim)
                if s > best:
                    best = s
    return best


# ---------------- one case
def check_case(twists, res, verbose=False):
    """res: dict with 'degrees', 'max_subsheaf_slope', 'kahler_einstein'."""
    k = len(twists)
    n = k + 2
    a, b, total, mu, F = prism_degrees(twists)
    intF = int_hex(F)
    # (a) degrees
    deg = [Fr(d) for d in res["degrees"]]
    ok_a = all(deg[i] == a[R6[i]] for i in range(6))
    ok_b = all(deg[6 + 2 * i] == b[i] and deg[6 + 2 * i + 1] == b[i]
               for i in range(k))
    ok_tot = sum(deg) == total
    # (b) reduced max slope vs harness
    ms = res["max_subsheaf_slope"]
    harness_max = Fr(ms[0], ms[1]) if isinstance(ms, list) else Fr(ms)
    red = reduced_max_slope(twists, a, b)
    ok_max = red == harness_max
    # (c) divergence identity: sum(a)/2 - mu == (n-1)!/2 * Phi
    Phi = Fr(0)
    for i, t in enumerate(twists):
        Fi = {(0, 0): Fr(1)}
        for j, s in enumerate(twists):
            if j != i:
                Fi = pmul(Fi, two_plus_ell(s))
        Phi += int_hex(pmul(Fi, ell_poly(t))) if ell_poly(t) else Fr(0)
    ok_div = (Fr(sum(a.values()), 2) - mu) == Fr(factorial(n - 1), 2) * Phi
    # (d) Full barycenter formula, including the KE moment criterion.
    m1 = int_hex(pmul(F, {(1, 0): Fr(1)}))
    m2 = int_hex(pmul(F, {(0, 1): Fr(1)}))
    fibre_bary = [m1 / intF, m2 / intF]
    expected_bary = fibre_bary + [
        -Fr(t[0] * m1 + t[1] * m2, 2 * intF) for t in twists
    ]
    stored_bary = [Fr(num, den) for num, den in res["barycenter"]]
    ok_bary = stored_bary == expected_bary and \
        ((m1 == 0 and m2 == 0) == bool(res["kahler_einstein"]))
    ok = ok_a and ok_b and ok_tot and ok_max and ok_div and ok_bary
    tag = "".join("." if x else "X" for x in
                  (ok_a, ok_b, ok_tot, ok_max, ok_div, ok_bary))
    if verbose or not ok:
        print(f"  {'OK ' if ok else 'FAIL'} [{tag}] "
              f"k={k} twists={[KEY[t] for t in twists]} "
              f"maxslope={harness_max} Phi={Phi}", flush=True)
    return ok


def main():
    # Keep the exact polynomial/integration utilities importable without
    # NumPy/SciPy; only this full-dimensional cross-check needs them.
    from root_twist_fast import CASES
    from toric_fast import analyze_fast

    random.seed(7)
    nz = [v for v in R6]
    allt = nz + [(0, 0)]
    total = passed = 0

    print("== canonical cases (from output/root_twist.json) ==", flush=True)
    with open(os.path.join(HERE, "..", "output", "root_twist.json")) as f:
        stored = {r["label"]: r for r in json.load(f)}
    twmap = {label: [VEC[w] for w in tw] for label, tw in CASES}
    for label, res in stored.items():
        total += 1
        passed += check_case(twmap[label], res, verbose=True)

    print("== all k=1 and k=2 tuples over roots u {0} ==", flush=True)
    todo = [[t] for t in allt] + [list(p) for p in iproduct(allt, allt)]
    for tw in todo:
        res = analyze_fast(root_twist_rays([KEY[t] for t in tw]))
        total += 1
        passed += check_case(tw, res)

    print("== random k=3 (24) and k=4 (10) tuples ==", flush=True)
    sample = [ [random.choice(allt) for _ in range(3)] for _ in range(24) ] \
           + [ [random.choice(allt) for _ in range(4)] for _ in range(10) ]
    for tw in sample:
        res = analyze_fast(root_twist_rays([KEY[t] for t in tw]))
        total += 1
        passed += check_case(tw, res)

    assert total == 98, total
    print(f"\nRESULT: {passed}/{total} cases pass ALL six checks "
          f"(degrees a, degrees b, total, reduced-max, divergence, barycenter/KE)",
          flush=True)
    if passed != total:
        raise SystemExit(1)
    print("PRISM LEMMA + REDUCTION + IDENTITIES: VALIDATED", flush=True)


if __name__ == "__main__":
    main()
