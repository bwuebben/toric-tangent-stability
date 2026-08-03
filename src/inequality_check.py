#!/usr/bin/env python3
"""Exact verification of the NEW inequality proofs (open-inequality session).

Verifies, in exact rational arithmetic, every step of the following claims
before they enter the paper:

CLAIM 1 (A' for negation-symmetric multisets).  If the twist multiset is
  invariant under t -> -t (so F is even), then for every twist t,
      int_Q F * l_t/(2+l_t)  =  - int_Q F * l_t^2/(4 - l_t^2)  <  0.

CLAIM 2 (A' for pure rotation-triples).  For F = H^m, H = prod over the
  triple (2+l_1)(2+l_2)(2+l_3), l_1+l_2+l_3 = 0:
      3 * int_Q l_1 (2+l_2)(2+l_3) H^{m-1}
        = int_Q (-2 p_2 + 3 e_3) H^{m-1}  <  0,
  via the pointwise domination 3|e_3| <= (3/2) p_2 and |odd part| <= even
  part of H^{m-1} = (A+B)^{m-1}, A = 8 - p_2, B = e_3.

CLAIM 3 (line inequality for X((a,-a)^r), exact margin).  With
  q = 4 - l_a^2, F = q^r, slicing Q over l_a = s (slice length 2-|s|):
      I_r = int_Q q^r = 4 K_r - 2 L_r,
      J_r' = int_Q q^{r-1} l_a^2 = 2[2(4K_{r-1}-K_r) - (4L_{r-1}-L_r)],
      K_r = int_0^1 (4-s^2)^r ds,   (2r+1) K_r = 3^r + 8r K_{r-1},
      L_r = int_0^1 s(4-s^2)^r ds = (4^{r+1}-3^{r+1})/(2(r+1)),
  and the critical inequality  2*3^r + 2r*J' < I_r  holds with margin
      I_r - 2*3^r - 2r*J'  =  (4^{r+1} - 3^{r+1})/(r+1)  >  0.

CLAIM 4 (full stability of the families, via the validated reduction).
  Using prism degrees + the reduced subspace list (both validated in
  prism_check.py against the harness), X((a,-a)^r) is stable for r = 1..6
  (dims 4..14) and X(triple^m) for m = 1..4 (dims 5..14) -- far beyond the
  complete census carried out in the paper for the larger cases.

CLAIM 5 is the line-share identity used in the reduction.  CLAIM 6 checks
every displayed integral identity in the proof for the unbalanced family
X((a,-a)^s,a,b), and verifies stability and non-KE exactly for s=0..20.
CLAIM 7 checks the six-sector identity and the extreme-exponent moment
ordering used in the zero-sum Kahler--Einstein classification.

Also reports the (B)-margins for the triple family and re-verifies KE
(moment = 0) for the symmetric families.
"""

import os
import sys
import json
from fractions import Fraction as Fr
from itertools import permutations
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from prism_check import (pmul, two_plus_ell, ell_poly, int_hex, int_edge,
                         prism_degrees, reduced_max_slope, R6)   # noqa: E402
from root_twist import root_twist_rays                           # noqa: E402
from toric_stability import det                                  # noqa: E402

A, B_, AB = (1, 0), (0, 1), (1, 1)
NEG = lambda t: (-t[0], -t[1])
TRIPLE = [A, B_, NEG(AB)]      # a, b, -(a+b): l_1 + l_2 + l_3 = 0


def prod_poly(multiset):
    F = {(0, 0): Fr(1)}
    for t in multiset:
        F = pmul(F, two_plus_ell(t))
    return F


def singleton(multiset, i):
    """int_Q F * l_{t_i}/(2+l_{t_i}) = int_Q l_{t_i} * prod_{j!=i}(2+l_{t_j})."""
    Fhat = prod_poly([t for j, t in enumerate(multiset) if j != i])
    return int_hex(pmul(ell_poly(multiset[i]), Fhat))


def check_claim1():
    print("== CLAIM 1: A' for negation-symmetric multisets ==")
    tests = [[A, NEG(A)], [A, NEG(A), B_, NEG(B_)],
             [A, NEG(A), AB, NEG(AB)],
             [A, A, NEG(A), NEG(A), B_, NEG(B_)],
             [A, NEG(A), B_, NEG(B_), AB, NEG(AB)]]
    ok = True
    for ms in tests:
        F = prod_poly(ms)
        for i, t in enumerate(ms):
            v = singleton(ms, i)
            # closed form: -int_Q F l^2/(4-l^2) = -int_Q (F/(2+l)/(2-l)) l^2
            Fhat2 = prod_poly([s for j, s in enumerate(ms)
                               if j != i and s != NEG(t)]
                              + [s for j, s in enumerate(ms)][:0])
            # build F/( (2+l_t)(2-l_t) ): remove one copy of t and one of -t
            rest = list(ms)
            rest.remove(t)
            rest.remove(NEG(t))
            l2 = pmul(ell_poly(t), ell_poly(t))
            closed = -int_hex(pmul(l2, prod_poly(rest)))
            ok &= (v == closed) and (v < 0)
        print(f"  ms={ms}: all singletons < 0 and == closed form: {ok}")
    print("  CLAIM 1:", "VERIFIED" if ok else "FAILED")
    return ok


def check_claim2(mmax=8):
    print("== CLAIM 2: A' for pure triples H^m ==")
    l1, l2, l3 = (ell_poly(t) for t in TRIPLE)
    H = prod_poly(TRIPLE)
    p2 = {}
    for l in (l1, l2, l3):
        for k, v in pmul(l, l).items():
            p2[k] = p2.get(k, Fr(0)) + v
    e3 = pmul(pmul(l1, l2), l3)
    ok = True
    Hm1 = {(0, 0): Fr(1)}          # H^{m-1}
    for m in range(1, mmax + 1):
        ms = TRIPLE * m
        v = singleton(ms, 0)       # t = gamma_1
        # cyclic identity: 3v == int (-2 p2 + 3 e3) H^{m-1}
        integrand = {}
        for k, c in p2.items():
            integrand[k] = integrand.get(k, Fr(0)) - 2 * c
        for k, c in e3.items():
            integrand[k] = integrand.get(k, Fr(0)) + 3 * c
        rhs = int_hex(pmul(integrand, Hm1))
        ok &= (3 * v == rhs) and (v < 0)
        print(f"  m={m}: singleton={v}  (3x == cyclic form: {3*v==rhs}, <0: {v<0})")
        Hm1 = pmul(Hm1, H)
    print("  CLAIM 2:", "VERIFIED" if ok else "FAILED")
    return ok


def check_claim3(rmax=40, cross2d=6):
    print("== CLAIM 3: line inequality for X((a,-a)^r), exact margin ==")
    K = [Fr(1)]                     # K_0 = 1
    L = [Fr(1, 2)]                  # L_0 = 1/2
    ok = True
    for r in range(1, rmax + 1):
        K.append(Fr(3 ** r + 8 * r * K[r - 1], 2 * r + 1))
        L.append(Fr(4 ** (r + 1) - 3 ** (r + 1), 2 * (r + 1)))
        I = 4 * K[r] - 2 * L[r]
        J = 2 * (2 * (4 * K[r - 1] - K[r]) - (4 * L[r - 1] - L[r]))
        margin = I - 2 * Fr(3 ** r) - 2 * r * J
        pred = Fr(4 ** (r + 1) - 3 ** (r + 1), r + 1)
        ok &= (margin == pred) and (margin > 0)
        # other lines: 2 K_r < I_r  <=>  L_r < K_r
        ok &= L[r] < K[r]
        if r <= cross2d:
            # 2-D cross-checks against the hexagon integrator
            q = {(0, 0): Fr(4), (2, 0): Fr(-1)}
            qr = {(0, 0): Fr(1)}
            for _ in range(r):
                qr = pmul(qr, q)
            ok &= int_hex(qr) == I
            qr1 = {(0, 0): Fr(1)}
            for _ in range(r - 1):
                qr1 = pmul(qr1, q)
            ok &= int_hex(pmul(qr1, {(2, 0): Fr(1)})) == J
            # edge values: q = 3 on E_{+-a}, edge length 1
            ok &= int_edge(qr, (1, 0)) == Fr(3 ** r)
            ok &= int_edge(qr, (-1, 0)) == Fr(3 ** r)
    print(f"  r=1..{rmax}: margin == (4^(r+1)-3^(r+1))/(r+1) > 0: {ok} "
          f"(2-D cross-checked r<={cross2d})")
    print("  CLAIM 3:", "VERIFIED" if ok else "FAILED")
    return ok


def stability_via_reduction(multiset, label):
    """Full reduced-criterion stability check from prism data alone."""
    a, b, total, mu, F = prism_degrees(multiset)
    best = reduced_max_slope(multiset, a, b)
    n = len(multiset) + 2
    verdict = ("stable" if best < mu else
               "strictly semistable" if best == mu else "unstable")
    # KE moment
    m1 = int_hex(pmul(F, {(1, 0): Fr(1)}))
    m2 = int_hex(pmul(F, {(0, 1): Fr(1)}))
    ke = (m1 == 0 and m2 == 0)
    print(f"  {label}: dim {n}, rho {len(multiset)+4}, (-K)^n={total}  "
          f"mu={mu}  max={best}  -> {verdict}  KE={ke}")
    return verdict, ke


def check_claim4():
    print("== CLAIM 4: full stability of the families (reduced criterion) ==")
    ok = True
    for r in range(1, 7):
        v, ke = stability_via_reduction([A, NEG(A)] * r, f"X((a,-a)^{r})")
        ok &= (v == "stable") and ke
    for m in range(1, 5):
        v, ke = stability_via_reduction(TRIPLE * m, f"X(triple^{m})")
        ok &= (v == "stable") and ke
    print("  CLAIM 4:", "VERIFIED" if ok else "FAILED")
    return ok


def triple_line_margins(mmax=8):
    print("== (B)-margins for X(triple^m) (exact symmetry check) ==")
    H = prod_poly(TRIPLE)
    Hm = dict(H)
    for m in range(1, mmax + 1):
        F = Hm
        intF = int_hex(F)
        worst = None
        for rho in [(1, 0), (0, 1), (-1, -1)]:
            edge = int_edge(F, rho) + int_edge(F, NEG(rho))
            corr = -m * singleton(TRIPLE * m,
                                  (TRIPLE * m).index(rho))  # -int F f_t > 0
            margin = intF - edge - corr
            if worst is None or margin < worst[0]:
                worst = (margin, rho)
        print(f"  m={m}: worst line-margin {worst[0]} "
              f"(= {float(worst[0]/intF):.4f} of int F)  positive: {worst[0]>0}")
        Hm = pmul(Hm, H)


def check_sum_identity():
    """CLAIM 5: given (a) [all singletons < 0], the three line-shares
    S_rho = int_{E_rho}F + int_{E_-rho}F + sum_{t_i=+-rho} (-varsigma_i)
    sum to EXACTLY 2 int_Q F  (divergence identity + every twist lies on
    exactly one line).  Hence (b) <=> strict triangle inequality."""
    print("== CLAIM 5: sum of line shares == 2 int_Q F ==")
    LINES = {(1, 0): [A, NEG(A)], (0, 1): [B_, NEG(B_)],
             (1, 1): [AB, NEG(AB)]}
    tests = [("dim7 counterexample", [A, B_, NEG(AB), A, NEG(A)]),
             ("pairs 2 lines", [A, NEG(A), B_, NEG(B_)]),
             ("triple", list(TRIPLE)),
             ("triple + pair (mixed)", TRIPLE + [B_, NEG(B_)]),
             ("triple^2 + pair", TRIPLE * 2 + [A, NEG(A)])]
    ok = True
    for label, ms in tests:
        F = prod_poly(ms)
        intF = int_hex(F)
        sings = [singleton(ms, i) for i in range(len(ms))]
        ok &= all(value < 0 for value in sings)
        assert all(s < 0 for s in sings), (label, sings)
        shares = {}
        for rho, roots in LINES.items():
            s = int_edge(F, rho) + int_edge(F, NEG(rho))
            s += sum(-sings[i] for i, t in enumerate(ms) if t in roots)
            shares[rho] = s
        tot = sum(shares.values())
        tri = all(2 * v < tot for v in shares.values())
        ok &= (tot == 2 * intF)
        print(f"  {label}: sum(shares)==2*intF: {tot == 2*intF}  "
              f"strict-triangle(=stable-side (b)): {tri}  "
              f"shares/intF={[str(v/intF) for v in shares.values()]}")
    print("  CLAIM 5:", "VERIFIED" if ok else "FAILED")
    return ok


def _moment_q(s, j):
    """M_j = int_0^1 u^j (4-u^2)^s du, exactly."""
    return sum(Fr(comb(s, h) * 4 ** (s - h) * (-1) ** h,
                  j + 2 * h + 1)
               for h in range(s + 1))


def check_claim6(smax=20):
    """Displayed formulas for Y_s = X((a,-a)^s,a,b)."""
    print("== CLAIM 6: unbalanced stable non-KE family ==")
    ok = True
    with open(os.path.join(os.path.dirname(HERE), "data",
                           "smooth_toric_fano_4d.json")) as f:
        census_rays = {p["id"]: p["vertices"] for p in
                       json.load(f)["polytopes"]}["F.4D.0061"]
    change = [[0, 0, -1, 0], [-1, 0, 1, 0],
              [0, 0, 0, 1], [0, 1, 0, 0]]
    mapped = {
        tuple(sum(v[i] * change[i][j] for i in range(4))
              for j in range(4))
        for v in root_twist_rays(["a", "b"])
    }
    identification_ok = abs(det(change)) == 1 and \
        mapped == {tuple(v) for v in census_rays}
    ok &= identification_ok
    print(f"  Y_0 == F.4D.0061 by an explicit GL(4,Z) map: "
          f"{identification_ok}")
    for s in range(smax + 1):
        ms = [A, NEG(A)] * s + [A, B_]
        F = prod_poly(ms)
        intF = int_hex(F)
        M = [_moment_q(s, j) for j in range(6)]
        ok &= intF == 16 * M[0] - 8 * M[1] - 2 * M[2] + M[3]

        phi_a = singleton(ms, 2 * s)
        phi_b = singleton(ms, 2 * s + 1)
        phi_ab = -2 * M[2] + M[3]
        ok &= phi_a == phi_b == phi_ab < 0
        if s:
            phi_minus_a = singleton(ms, 1)
            predicted = -(24 * _moment_q(s - 1, 2)
                          - 12 * _moment_q(s - 1, 3)
                          - 2 * _moment_q(s - 1, 4)
                          + _moment_q(s - 1, 5))
            ok &= phi_minus_a == predicted < 0

        sings = [singleton(ms, i) for i in range(len(ms))]
        shares = {}
        for rho in (A, B_, AB):
            share = int_edge(F, rho) + int_edge(F, NEG(rho))
            share += sum(-sings[i] for i, t in enumerate(ms)
                         if t in (rho, NEG(rho)))
            shares[rho] = share
        ok &= shares[B_] == 8 * M[0] - 2 * M[1] + 2 * M[2] - M[3]
        ok &= shares[AB] == 8 * M[0] + 2 * M[1] - 2 * M[2]
        ok &= all(share < intF for share in shares.values())

        moment_x = int_hex(pmul(F, {(1, 0): Fr(1)}))
        predicted_moment = 2 * (2 * M[2] - M[3])
        ok &= moment_x == predicted_moment > 0
        # The singleton and share checks are precisely Proposition 7.5.
        # Cross-check its brute-force reduced-subspace implementation only
        # for the first few members; that implementation grows exponentially
        # with the number of twists.
        if s <= 3:
            verdict, ke = stability_via_reduction(ms, f"Y_{s}")
            ok &= verdict == "stable" and not ke
    print(f"  s=0..{smax}: all identities, stability, and non-KE: {ok}")
    print("  CLAIM 6:", "VERIFIED" if ok else "FAILED")
    return ok


def check_claim7(dmax=3, qmax=4):
    """KE necessity: a maximal pair exponent has smaller moment than a minimal one."""
    print("== CLAIM 7: zero-sum KE criterion and orbit ordering ==")
    ok = True

    # Verify the six-permutation identity independently on rational sectors.
    sectors = [(Fr(1, 7), Fr(2, 7)),
               (Fr(1, 5), Fr(3, 10)),
               (Fr(2, 9), Fr(1, 3))]
    for a, b in sectors:
        c = a + b
        vals = (4 - a * a, 4 - b * b, 4 - c * c)
        scores = (-a, -b, c)
        for q1 in range(qmax + 1):
            for q2 in range(qmax + 1):
                for q3 in range(qmax + 1):
                    qs = (q1, q2, q3)
                    if q1 == q2 == q3:
                        continue
                    i = max(range(3), key=qs.__getitem__)
                    j = min(range(3), key=qs.__getitem__)
                    k = 3 - i - j

                    def orbit_sum(coord):
                        return sum(
                            scores[p[coord]]
                            * vals[p[0]] ** qs[0]
                            * vals[p[1]] ** qs[1]
                            * vals[p[2]] ** qs[2]
                            for p in permutations(range(3))
                        )

                    q = qs[j]
                    r = qs[i] - q
                    h = qs[k] - q
                    aa, bb, cc = vals
                    rhs = (aa * bb * cc) ** q * (
                        (b - a) * (aa ** r - bb ** r) * cc ** h
                        - (2 * a + b) * (aa ** r - cc ** r) * bb ** h
                        - (a + 2 * b) * (bb ** r - cc ** r) * aa ** h
                    )
                    ok &= orbit_sum(i) - orbit_sum(j) == rhs < 0

    # Direct exact hexagon integration, independent of the sector algebra.
    for d in range(1, dmax + 1):
        for q1 in range(qmax + 1):
            for q2 in range(qmax + 1):
                for q3 in range(qmax + 1):
                    qs = (q1, q2, q3)
                    ms = (TRIPLE * d + [A, NEG(A)] * q1
                          + [B_, NEG(B_)] * q2
                          + [NEG(AB), AB] * q3)
                    F = prod_poly(ms)
                    moments = tuple(int_hex(pmul(F, ell_poly(t)))
                                    for t in TRIPLE)
                    ok &= sum(moments) == 0
                    if q1 == q2 == q3:
                        ok &= moments == (0, 0, 0)
                    else:
                        i = max(range(3), key=qs.__getitem__)
                        j = min(range(3), key=qs.__getitem__)
                        ok &= moments[i] < moments[j]

    # The d=0 density is centrally even for arbitrary pair exponents.
    for q1 in range(qmax + 1):
        for q2 in range(qmax + 1):
            for q3 in range(qmax + 1):
                ms = ([A, NEG(A)] * q1 + [B_, NEG(B_)] * q2
                      + [NEG(AB), AB] * q3)
                F = prod_poly(ms)
                ok &= all(int_hex(pmul(F, ell_poly(t))) == 0
                          for t in TRIPLE)

    print(f"  orbit identity + exact grid d=1..{dmax}, "
          f"q_i=0..{qmax}: {ok}")
    print("  CLAIM 7:", "VERIFIED" if ok else "FAILED")
    return ok


def main():
    r1 = check_claim1()
    r2 = check_claim2()
    r3 = check_claim3()
    r4 = check_claim4()
    r5 = check_sum_identity()
    r6 = check_claim6()
    r7 = check_claim7()
    triple_line_margins()
    ok = all((r1, r2, r3, r4, r5, r6, r7))
    print("\nALL PROOF CLAIMS:", "VERIFIED" if ok else "SOMETHING FAILED")
    if not ok:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
