#!/usr/bin/env python3
"""Enumerate ALL zero-sum multisets of nonzero roots up to the hexagon's
dihedral symmetry, for k = 2..6 twists (dimensions 4 through 8), and verify
on each class: (i) T_X is (-K)-stable, and (ii) X is Kahler-Einstein exactly
when the twist multiset is invariant under negation or under the order-three
rotation of the root hexagon.

This regenerates, using the exact rational prism formulas, the classification
behind the Kahler-Einstein criterion table of Section 7 of the paper.

Usage:  python3 roottwist_classes.py          # writes output/root_twist_classes.json
"""

import json
import os
from itertools import combinations_with_replacement

from prism_check import (ell_poly, int_hex, pmul, prism_degrees,
                         reduced_max_slope)
from root_twist import ROOTS

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "output")
JSON_PATH = os.path.join(OUT, "root_twist_classes.json")

NONZERO = [k for k in ROOTS if k != "0"]          # the six roots
VEC = {k: ROOTS[k] for k in NONZERO}
KEY_OF = {v: k for k, v in VEC.items()}           # vector -> key


def hexagon_symmetries():
    """All lattice automorphisms preserving the six-root set (dihedral, order 12)."""
    roots = list(VEC.values())
    maps = []
    for ia in roots:
        for ib in roots:
            # candidate M sends a=(1,0)->ia, b=(0,1)->ib; M is det +-1 with
            # M(a+b) = ia+ib also a root, and M must permute the root set.
            s = (ia[0] + ib[0], ia[1] + ib[1])
            if s not in KEY_OF:
                continue
            img = {(1, 0): ia, (0, 1): ib}
            def apply(v, img=img):
                return (v[0] * img[(1, 0)][0] + v[1] * img[(0, 1)][0],
                        v[0] * img[(1, 0)][1] + v[1] * img[(0, 1)][1])
            if {apply(r) for r in roots} == set(roots):
                maps.append(apply)
    # deduplicate by action on the roots
    seen, out = set(), []
    for m in maps:
        sig = tuple(m(r) for r in sorted(VEC.values()))
        if sig not in seen:
            seen.add(sig)
            out.append(m)
    return out


SYMS = hexagon_symmetries()
assert len(SYMS) == 12, f"expected dihedral group of order 12, got {len(SYMS)}"

NEG = lambda v: (-v[0], -v[1])
TAU = None
for m in SYMS:  # order-three rotation: a -> b -> -(a+b)
    if m(VEC["a"]) == VEC["b"] and m(VEC["b"]) == VEC["-(a+b)"]:
        TAU = m
assert TAU is not None


def canon(multiset_keys):
    """Canonical representative of the multiset class under the dihedral group."""
    best = None
    for m in SYMS:
        img = tuple(sorted(KEY_OF[m(VEC[k])] for k in multiset_keys))
        if best is None or img < best:
            best = img
    return best


def invariant(multiset_keys, g):
    img = sorted(KEY_OF[g(VEC[k])] for k in multiset_keys)
    return img == sorted(multiset_keys)


def main():
    classes = {}
    for k in range(2, 7):
        for combo in combinations_with_replacement(NONZERO, k):
            sx = sum(VEC[t][0] for t in combo)
            sy = sum(VEC[t][1] for t in combo)
            if (sx, sy) != (0, 0):
                continue
            classes.setdefault(canon(combo), combo)
    by_k = {}
    for rep in classes.values():
        by_k.setdefault(len(rep), []).append(rep)
    counts = {k: len(v) for k, v in sorted(by_k.items())}
    print(f"zero-sum classes by k: {counts}  (total {len(classes)})")
    assert counts == {2: 1, 3: 1, 4: 2, 5: 1, 6: 4}, counts

    results = []
    for rep in sorted(classes.values(), key=lambda r: (len(r), r)):
        n = 2 + len(rep)
        label = f"n={n} ({','.join(rep)})"
        neg_inv = invariant(rep, NEG)
        tau_inv = invariant(rep, TAU)
        twists = [VEC[t] for t in rep]
        a, b, total, mu, F = prism_degrees(twists)
        max_slope = reduced_max_slope(twists, a, b)
        moment = [int_hex(pmul(F, ell_poly(v))) for v in ((1, 0), (0, 1))]
        verdict = ("stable" if max_slope < mu else
                   "strictly semistable" if max_slope == mu else "unstable")
        ke = all(v == 0 for v in moment)
        row = {
            "n": n, "twists": list(rep),
            "negation_invariant": neg_inv, "rotation_invariant": tau_inv,
            "verdict": verdict, "kahler_einstein": ke,
            "anticanonical_degree": int(total),
            "mu_TX": [mu.numerator, mu.denominator],
            "max_subsheaf_slope": [max_slope.numerator, max_slope.denominator],
        }
        results.append(row)
        print(f"CLASS {label}: {verdict}  KE={ke}  neg={neg_inv} "
              f"tau={tau_inv}  (-K)^n={int(total)}",
              flush=True)
        assert verdict == "stable", (label, verdict)
        assert ke == (neg_inv or tau_inv), (label, ke, neg_inv, tau_inv)

    with open(JSON_PATH, "w") as f:
        json.dump(results, f, indent=1)
    ok = sum(1 for r in results if r["kahler_einstein"])
    print(f"\nRESULT: {len(results)}/9 zero-sum classes stable; "
          f"KE in exactly the {ok} classes with a symmetry; criterion 0/{len(results)} mismatches")
    print(f"wrote {JSON_PATH}")


if __name__ == "__main__":
    main()
