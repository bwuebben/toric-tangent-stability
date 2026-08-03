#!/usr/bin/env python3
"""Polystability refinement for the strictly-semistable tangent bundles.

An equivariant splitting T_X = F_{V_1} + ... + F_{V_m} exists iff N_C = (+) V_i
with every ray line C.u_rho contained in some V_i (Klyachko filtrations are
two-step, so a direct-sum decomposition of the fiber compatible with all
filtrations is exactly a partition of the rays into groups whose spans are
linearly independent and exhaust N).

T_X is POLYSTABLE iff it is semistable and such a decomposition exists in
which every piece F_{V_i} is stable of slope mu(T_X):
  * slope(F_{V_i}) = sum_{rho in V_i} deg(D_rho) / dim V_i = mu, and
  * every proper nonzero ray-spanned subspace U of V_i has slope(U) < mu.

By Kobayashi--Hitchin, polystability is equivalent to the existence of a
Hermitian--Einstein metric on T_X with respect to a Kahler form representing
c_1(X).  So this module separates the strictly-semistable class into HE
(polystable) and non-HE cases.

Usage: python3 polystability.py 3   (reads output/sweep_3d.json + data JSON)
"""

import json
import os
import sys
from fractions import Fraction as Fr
from itertools import combinations

from toric_stability import rref, in_span

HERE = os.path.dirname(os.path.abspath(__file__))


def ray_closed_subspaces(rays):
    """All subspaces spanned by rays, as (member frozenset, dim, rref)."""
    n = len(rays[0])
    r = len(rays)
    seen = {}
    for k in range(1, n + 1):
        for comb in combinations(range(r), k):
            rows = [rays[i] for i in comb]
            rk, canon = rref(rows)
            if rk < k or canon in seen:
                continue
            members = frozenset(i for i in range(r) if in_span(canon, rays[i]))
            seen[canon] = (members, k)
    return list(seen.values())


def stable_pieces_of_slope_mu(rays, deg, mu):
    """Candidate direct summands: ray-closed subspaces V with slope(V) = mu
    and every proper ray-spanned subspace strictly below mu."""
    pieces = []
    for members, dim in ray_closed_subspaces(rays):
        if Fr(sum(deg[i] for i in members), dim) != mu:
            continue
        sub_rays = [rays[i] for i in sorted(members)]
        idx = sorted(members)
        stable = True
        for smem, sdim in ray_closed_subspaces(sub_rays):
            if sdim == dim:
                continue
            if Fr(sum(deg[idx[j]] for j in smem), sdim) >= mu:
                stable = False
                break
        if stable:
            pieces.append((members, dim))
    return pieces


def is_polystable(rays, deg, mu):
    """Exact-cover search: partition all rays into stable slope-mu pieces
    whose spans are independent and exhaust N.  Returns a decomposition
    (list of ray-index sets) or None."""
    n = len(rays[0])
    r = len(rays)
    pieces = stable_pieces_of_slope_mu(rays, deg, mu)

    def search(assigned, dim_used, chosen):
        if len(assigned) == r:
            return list(chosen) if dim_used == n else None
        pivot = min(set(range(r)) - assigned)
        for members, dim in pieces:
            if pivot not in members or members & assigned:
                continue
            if dim_used + dim > n:
                continue
            # independence: combined ray set must have rank dim_used + dim
            rows = [rays[i] for g in chosen for i in g] + \
                   [rays[i] for i in sorted(members)]
            if rref(rows)[0] != dim_used + dim:
                continue
            chosen.append(sorted(members))
            res = search(assigned | members, dim_used + dim, chosen)
            if res is not None:
                return res
            chosen.pop()
        return None

    return search(frozenset(), 0, [])


def ray_matroid_components(rays):
    """Connected components of the linear matroid represented by ``rays``.

    Fundamental circuits relative to any basis generate the matroid's
    connectivity relation.  The result is a canonical partition of the ray
    indices, independent of the greedy basis used to compute it.
    """
    r = len(rays)
    n = len(rays[0])

    basis = []
    for i, ray in enumerate(rays):
        if rref([rays[j] for j in basis] + [ray])[0] > len(basis):
            basis.append(i)
        if len(basis) == n:
            break
    assert len(basis) == n, "ray generators do not span the ambient lattice"

    parent = list(range(r))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    def union(i, j):
        i, j = find(i), find(j)
        if i != j:
            parent[j] = i

    basis_set = set(basis)
    for e in range(r):
        if e in basis_set:
            continue
        circuit = [e]
        for b in basis:
            other_basis = [rays[j] for j in basis if j != b]
            canon = () if not other_basis else rref(other_basis)[1]
            if not in_span(canon, rays[e]):
                circuit.append(b)
        assert len(circuit) >= 2
        for i in circuit[1:]:
            union(circuit[0], i)

    groups = {}
    for i in range(r):
        groups.setdefault(find(i), []).append(i)
    return sorted((sorted(group) for group in groups.values()),
                  key=lambda group: group[0])


def is_polystable_matroid(rays, deg, mu):
    """Canonical ray-matroid test for polystability.

    Returns the component partition when every component sheaf is stable of
    slope ``mu``, and ``None`` otherwise.
    """
    components = ray_matroid_components(rays)
    for component in components:
        component_rays = [rays[i] for i in component]
        dim = rref(component_rays)[0]
        if Fr(sum(deg[i] for i in component), dim) != mu:
            return None
        for members, subdim in ray_closed_subspaces(component_rays):
            if subdim == dim:
                continue
            slope = Fr(sum(deg[component[j]] for j in members), subdim)
            if slope >= mu:
                return None
    return components


def main(dim):
    with open(os.path.join(HERE, "..", "output", f"sweep_{dim}d.json")) as f:
        results = json.load(f)
    with open(os.path.join(HERE, "..", "data",
                           f"smooth_toric_fano_{dim}d.json")) as f:
        db = {p["id"]: p for p in json.load(f)["polytopes"]}

    out = []
    n_poly = n_not = 0
    for res in results:
        if res["verdict"] != "strictly semistable":
            continue
        rays = [tuple(v) for v in db[res["name"]]["vertices"]]
        deg = res["degrees"]
        mu = Fr(res["mu_TX"][0], res["mu_TX"][1])
        decomp = is_polystable(rays, deg, mu)
        matroid_decomp = is_polystable_matroid(rays, deg, mu)
        assert (decomp is not None) == (matroid_decomp is not None), \
            f"splitting tests disagree for {res['name']}"
        ok = decomp is not None
        n_poly += ok
        n_not += (not ok)
        out.append({"id": res["name"], "polystable": ok,
                    "decomposition_ray_indices": decomp,
                    "KE": res["kahler_einstein"]})
        tag = "POLYSTABLE (HE exists)" if ok else "not polystable (no HE)"
        ke = "KE " if res["kahler_einstein"] else "    "
        print(f"#{res['name']:<14} {ke} {tag}"
              + (f"   ray partition: {decomp}" if ok else ""), flush=True)

    path = os.path.join(HERE, "..", "output", f"polystability_{dim}d.json")
    with open(path, "w") as f:
        json.dump(out, f, indent=1)
    print(f"\nstrictly semistable in dim {dim}: {len(out)} "
          f"-> polystable {n_poly}, not polystable {n_not}")
    print(f"wrote {os.path.normpath(path)}")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
