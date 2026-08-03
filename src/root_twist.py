#!/usr/bin/env python3
"""The root-twisted dP3-fibration family.

X(t_1, ..., t_k): a toric dP3-fibration over (P^1)^k, dimension n = 2 + k.
Fan rays in Z^{2+k}:
  * the A2 hexagon in the fiber plane:  +-a, +-b, +-(a+b),  a = e1, b = e2;
  * for each base factor i = 1..k:  f_i + t_i  and  -f_i,  f_i = e_{2+i},
    with twist t_i a vector in the fiber plane (typically one of the six
    roots {+-a, +-b, +-(a+b)} or 0).

t_i = 0 for all i gives dP3 x (P^1)^k (polystable, KE).  F.5D.0611 is
(up to lattice iso) k = 3 with twists (a, b, -(a+b)).

Theorems under regression: nonzero root twists summing to zero give a stable
toric Fano; among them, the variety is Kahler--Einstein exactly for a
negation- or rotation-invariant twist multiset.  Zero sum is not necessary
for stability: X((a,-a)^s,a,b) is a proved stable non-KE family.  Cases of
dimension 7 and 8 lie beyond the complete census in the paper.

Usage: python3 root_twist.py
"""

import sys
from fractions import Fraction as Fr

from toric_stability import analyze, fmt_frac

A = (1, 0)
B = (0, 1)
AB = (1, 1)


def neg(v):
    return tuple(-x for x in v)


ROOTS = {"a": A, "b": B, "a+b": AB, "-a": neg(A), "-b": neg(B),
         "-(a+b)": neg(AB), "0": (0, 0)}


def root_twist_rays(twists):
    """Rays of X(t_1..t_k) in Z^{2+k}; twists = list of keys into ROOTS."""
    k = len(twists)
    n = 2 + k

    def emb(fiber2, base=None, sign=1):
        v = [0] * n
        v[0], v[1] = fiber2
        if base is not None:
            v[2 + base] = sign
        return tuple(v)

    rays = [emb(A), emb(neg(A)), emb(B), emb(neg(B)), emb(AB), emb(neg(AB))]
    for i, t in enumerate(twists):
        rays.append(emb(ROOTS[t], base=i, sign=1))     # f_i + t_i
        rays.append(emb((0, 0), base=i, sign=-1))      # -f_i
    return rays


def main():
    cases = [
        ("k=3  (a, b, -(a+b))     [should be F.5D.0611]",
         ["a", "b", "-(a+b)"]),
        ("k=3  (0, 0, 0)          [product dP3 x (P1)^3]",
         ["0", "0", "0"]),
        ("k=4  (a, -a, b, -b)     [zero-sum pairs, dim 6]",
         ["a", "-a", "b", "-b"]),
        ("k=4  (a, b, -(a+b), 0)  [one untwisted factor, dim 6]",
         ["a", "b", "-(a+b)", "0"]),
        ("k=4  (a, b, -(a+b), a)  [sum = a != 0, dim 6]",
         ["a", "b", "-(a+b)", "a"]),
        ("k=5  (a, b, -(a+b), a, -a)   [zero-sum, dim 7]",
         ["a", "b", "-(a+b)", "a", "-a"]),
        ("k=6  (a, b, -(a+b)) x2       [zero-sum, dim 8]",
         ["a", "b", "-(a+b)", "a", "b", "-(a+b)"]),
    ]
    for name, twists in cases:
        rays = root_twist_rays(twists)
        try:
            res = analyze(rays, name)
        except AssertionError as e:
            print(f"{name:<44} NOT smooth/reflexive: {e}")
            continue
        print(f"{name:<44} dim {res['dim']}  (-K)^n = "
              f"{res['anticanonical_degree']:>5}  mu(T) = "
              f"{fmt_frac(res['mu_TX']):>8}  max mu(V) = "
              f"{fmt_frac(res['max_subsheaf_slope']):>8}  "
              f"{res['verdict']:<20} KE: "
              f"{'yes' if res['kahler_einstein'] else 'no'}", flush=True)


if __name__ == "__main__":
    main()
