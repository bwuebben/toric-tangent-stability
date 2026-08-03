#!/usr/bin/env python3
"""Faster hybrid stability analysis for the high-dimensional root-twist
cases where the naive dual-vertex enumeration in toric_stability.analyze is the
bottleneck (it exact-solves all C(r,n) ray n-subsets).

Only the dual-polytope VERTEX ENUMERATION is replaced: a numpy float pre-filter
(det + feasibility with a loose tolerance) discards the ~99% of n-subsets that
are singular or infeasible, and every surviving candidate is confirmed in EXACT
integer arithmetic.  Exact confirmation rules out false-positive vertices,
but floating-point candidate generation and QHull triangulation are not, by
themselves, certificates against false negatives or an invalid triangulation.
Accordingly this module is an accelerated cross-check, not the source of
record for the paper's high-dimensional claims.  For the reported
inputs, an exhaustive exact basis audit found no omitted vertices, and the
independent exact two-dimensional prism formulas in prism_check.py reproduce
all degrees, slopes, and barycenters.

Each survivor is verified exactly: the rounded integer point must satisfy its
n tight equations <m,u_i> = -1 and all halfspaces <m,u> >= -1.

Everything downstream (facet volumes, barycenter, the subspace slope loop) reuses
the exact primitives of toric_stability unchanged.  `_selftest` cross-checks
analyze_fast against toric_stability.analyze on small cases.
"""

import math
from fractions import Fraction as Fr
from itertools import combinations

import numpy as np
from scipy.spatial import Delaunay

import toric_stability as ts
from toric_stability import (dot, rref, in_span, det as _det,
                             _to_lattice_coords,
                             barycenter_and_volume_from_facets)


def facet_volume_and_triangulation_fast(u, facet_verts):
    """Evaluate a QHull-proposed triangulation with exact determinants.

    The returned sum is exact conditional on the proposed simplices forming a
    valid triangulation; callers needing a certificate should use the exact
    prism formulas or the reference implementation.
    """
    coords, lift = _to_lattice_coords(u, facet_verts)
    d = len(coords[0])
    if d == 0:
        return 1, [(lift[coords[0]],)]
    if d == 1:
        lo, hi = min(coords), max(coords)
        return hi[0] - lo[0], [(lift[lo], lift[hi])]
    tri = Delaunay(np.array(coords, dtype=float))
    vol = 0
    lifted = []
    for simplex in tri.simplices:
        spts = [coords[i] for i in simplex]        # d+1 lattice points
        edges = [[spts[i][k] - spts[0][k] for k in range(d)]
                 for i in range(1, len(spts))]
        v = abs(int(_det(edges)))
        if v == 0:
            continue                               # coplanar sliver: 0 volume
        vol += v
        lifted.append(tuple(lift[p] for p in spts))
    return vol, lifted


def dual_polytope_vertices_fast(rays):
    """Same output as toric_stability.dual_polytope_vertices, float-accelerated."""
    n = len(rays[0])
    r = len(rays)
    R = np.array(rays, dtype=float)          # r x n
    bvec = -np.ones(n)
    verts = set()
    tried = set()
    for comb in combinations(range(r), n):
        A = R[list(comb)]
        det = np.linalg.det(A)
        if abs(det) < 0.5:                   # integer det: |det|>=1 or exactly 0
            continue
        x = np.linalg.solve(A, bvec)
        if not np.all(R @ x >= -1.0 - 1e-6): # loose: keep every true vertex
            continue
        key = tuple(int(v) for v in np.rint(x))
        if key in tried:
            continue
        tried.add(key)
        # exact confirmation
        if all(dot(key, u) == -1 for u in (rays[i] for i in comb)) \
           and all(dot(key, u) >= -1 for u in rays):
            verts.add(key)
    return sorted(verts)


def analyze_fast(rays, name="", log=None):
    """Hybrid accelerated verdict, with exact evaluation of retained data."""
    def emit(m):
        if log:
            log(m)
    n = len(rays[0])
    r = len(rays)
    emit(f"{name}: enumerating dual vertices (dim {n}, {r} rays)")
    dual_verts = dual_polytope_vertices_fast(rays)
    emit(f"{name}: {len(dual_verts)} dual vertices; computing facet volumes")
    deg = []
    facet_tris = []
    for u in rays:
        fverts = [m for m in dual_verts if dot(m, u) == -1]
        assert len(fverts) >= n, f"ray {u}: facet has too few vertices"
        v, tris = facet_volume_and_triangulation_fast(u, fverts)
        deg.append(v)
        facet_tris.append(tris)
    total = sum(deg)
    mu = Fr(total, n)
    bary, vol = barycenter_and_volume_from_facets(n, facet_tris)
    assert total == math.factorial(n) * vol, (name, total, vol)

    emit(f"{name}: (-K)^n={total} mu={mu}; scanning ray-spanned subspaces")
    seen = set()
    worst = (Fr(-1), None)
    top_canons = []
    equalities = []
    for k in range(1, n):
        cnt = 0
        for comb in combinations(range(r), k):
            rows = [rays[i] for i in comb]
            rk, canon = rref(rows)
            if rk < k or canon in seen:
                continue
            seen.add(canon)
            cnt += 1
            members = [i for i in range(r) if in_span(canon, rays[i])]
            slope = Fr(sum(deg[i] for i in members), k)
            desc = {"dim": k, "rays_in_V": [rays[i] for i in members],
                    "deg": sum(deg[i] for i in members),
                    "slope": [slope.numerator, slope.denominator]}
            if slope > worst[0]:
                worst = (slope, desc)
                top_canons = [canon]
            elif slope == worst[0]:
                top_canons.append(canon)
            if slope == mu:
                equalities.append(desc)
        emit(f"{name}:   level k={k}/{n-1} done ({cnt} distinct subspaces)")

    if worst[0] > mu:
        verdict = "unstable"
        hn_dim, hn_canon = rref(
            [row for canon in top_canons for row in canon])
        hn_members = [i for i in range(r) if in_span(hn_canon, rays[i])]
        hn_deg = sum(deg[i] for i in hn_members)
        hn_slope = Fr(hn_deg, hn_dim)
        assert hn_dim < n and hn_slope == worst[0]
        assert all(all(in_span(hn_canon, row) for row in canon)
                   for canon in top_canons)
        worst = (hn_slope,
                 {"dim": hn_dim,
                  "rays_in_V": [rays[i] for i in hn_members],
                  "deg": hn_deg,
                  "slope": [hn_slope.numerator, hn_slope.denominator]})
    elif worst[0] == mu:
        verdict = "strictly semistable"
    else:
        verdict = "stable"
    is_ke = all(c == 0 for c in bary)
    return {"name": name, "dim": n, "num_rays": r, "degrees": deg,
            "anticanonical_degree": total,
            "mu_TX": [mu.numerator, mu.denominator],
            "max_subsheaf_slope": worst[0] if worst[1] is None
            else worst[1]["slope"],
            "witness": worst[1], "equality_witnesses": equalities,
            "verdict": verdict, "barycenter": [[c.numerator, c.denominator]
                                               for c in bary],
            "kahler_einstein": is_ke}


def _selftest():
    """Cross-check analyze_fast == toric_stability.analyze on cases both can do
    quickly.  Higher-dimensional records are checked against the exact prism
    formulas by prism_check.py."""
    from root_twist import root_twist_rays
    cases = [("dim4-pair", ["a", "-a"]),
             ("dim4-unbalanced", ["a", "b"]),
             ("dim4-unstable", ["a", "a"]),
             ("dim5", ["a", "b", "-(a+b)"]),
             ("dim5-prod", ["0", "0", "0"]),
             ]
    for nm, tw in cases:
        rays = root_twist_rays(tw)
        a = ts.analyze(rays, nm)
        b = analyze_fast(rays, nm)
        for key in ("degrees", "anticanonical_degree", "mu_TX",
                    "max_subsheaf_slope", "verdict", "barycenter",
                    "kahler_einstein"):
            assert a[key] == b[key], (nm, key, a[key], b[key])
        # dual vertices identical
        assert (sorted(ts.dual_polytope_vertices(rays))
                == dual_polytope_vertices_fast(rays)), nm
        print(f"  OK {nm}: {b['verdict']} KE={b['kahler_einstein']} "
              f"(-K)^n={b['anticanonical_degree']}")
    print("selftest passed: analyze_fast matches analyze exactly")


if __name__ == "__main__":
    _selftest()
