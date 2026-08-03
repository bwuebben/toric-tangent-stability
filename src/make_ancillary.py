#!/usr/bin/env python3
"""Build the arXiv ancillary-file package (anc/) promised by Theorem 4.1.

Per dimension d = 3,4,5,6 emits anc/stability_<d>d.json.gz: a JSON array with
one record per smooth toric Fano d-fold:

  id            polyDB identifier (Polytopes.Lattice.SmoothReflexive)
  picard_rank
  rays          primitive ray generators of the fan (= vertices of the Fano
                polytope), in the fixed order used by all aligned fields
  degrees       anticanonical degrees deg_{-K} D_rho, aligned with rays
  anticanonical_degree      (-K)^d = sum(degrees)
  mu            slope of T_X, exact rational [num, den]
  max_slope     maximal slope over ray-spanned proper subspaces [num, den]
  witness       a subspace attaining max_slope: {dim, rays_in_V, deg, slope};
                for unstable records this is the first HN subspace
  verdict       stable | strictly semistable | unstable
  barycenter    exact barycenter of the moment polytope, [[num,den],...]
  kahler_einstein           (barycenter == 0)
  polystable    (strictly semistable only) does T_X split equivariantly into
                stable equal-slope summands
  splitting     (polystable only) the certifying partition of ray indices

Also emits anc/root_twist.json.gz (the Section-7 family, dims 4-8) plus
anc/root_twist_classes.json.gz (the nine zero-sum classes) and
anc/README.txt (schema + provenance + the consistency tallies).

Before writing each artifact, the script *asserts* its counts against the
numbers printed in the paper: classification counts (Thm 4.1), polystable
counts (Thm 5.4), and the KE column (Section 6 table).
"""

import gzip
import json
import os

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "output")
DATA = os.path.join(HERE, "..", "data")
ANC = os.path.join(HERE, "..", "paper", "anc")

# the numbers printed in the paper -- the package must reproduce them exactly
PAPER = {
    3: dict(total=18, stable=4, sss=7, unstable=7, poly=4, ke_stable=2),
    4: dict(total=124, stable=26, sss=24, unstable=74, poly=14, ke_stable=4),
    5: dict(total=866, stable=109, sss=82, unstable=675, poly=52, ke_stable=7),
    6: dict(total=7622, stable=636, sss=366, unstable=6620, poly=259,
            ke_stable=13),
}


def write_gzip_json(path, value):
    """Write deterministic gzip JSON (fixed timestamp, compact encoding)."""
    raw = json.dumps(value, separators=(",", ":")).encode("utf-8")
    with open(path, "wb") as f:
        f.write(gzip.compress(raw, compresslevel=9, mtime=0))


def build_dim(d):
    with open(os.path.join(OUT, f"sweep_{d}d.json")) as f:
        sweep = json.load(f)
    with open(os.path.join(DATA, f"smooth_toric_fano_{d}d.json")) as f:
        db = json.load(f)
    rays_by_id = {str(p["id"]): p["vertices"] for p in db["polytopes"]}
    with open(os.path.join(OUT, f"polystability_{d}d.json")) as f:
        poly = {r["id"]: r for r in json.load(f)}

    recs = []
    tally = {"stable": 0, "strictly semistable": 0, "unstable": 0}
    npoly = ke_stable = 0
    for r in sweep:
        vid = r["name"]
        rec = {
            "id": vid,
            "picard_rank": r["extra"]["picard_rank"],
            "rays": rays_by_id[vid],
            "degrees": r["degrees"],
            "anticanonical_degree": r["anticanonical_degree"],
            "mu": r["mu_TX"],
            "max_slope": r["max_subsheaf_slope"],
            "witness": r["witness"],
            "verdict": r["verdict"],
            "barycenter": r["barycenter"],
            "kahler_einstein": r["kahler_einstein"],
        }
        assert len(rec["rays"]) == len(rec["degrees"]) == r["num_rays"], vid
        assert sum(rec["degrees"]) == rec["anticanonical_degree"], vid
        tally[r["verdict"]] += 1
        if r["verdict"] == "stable" and r["kahler_einstein"]:
            ke_stable += 1
        if r["verdict"] == "strictly semistable":
            p = poly[vid]
            rec["polystable"] = bool(p["polystable"])
            if p["polystable"]:
                rec["splitting"] = p["decomposition_ray_indices"]
                npoly += 1
            assert bool(p["KE"]) == rec["kahler_einstein"], vid
        recs.append(rec)

    want = PAPER[d]
    assert len(recs) == want["total"], (d, len(recs))
    assert tally["stable"] == want["stable"], (d, tally)
    assert tally["strictly semistable"] == want["sss"], (d, tally)
    assert tally["unstable"] == want["unstable"], (d, tally)
    assert npoly == want["poly"], (d, npoly)
    assert ke_stable == want["ke_stable"], (d, ke_stable)

    path = os.path.join(ANC, f"stability_{d}d.json.gz")
    write_gzip_json(path, recs)
    print(f"dim {d}: {len(recs)} records "
          f"({tally['stable']}/{tally['strictly semistable']}/"
          f"{tally['unstable']}; polystable {npoly}; KE-stable {ke_stable}) "
          f"-> {os.path.basename(path)} "
          f"[{os.path.getsize(path)//1024} KB]")
    return tally, npoly


def build_root_twist():
    with open(os.path.join(OUT, "root_twist.json")) as f:
        rt = json.load(f)
    keep = ["label", "verdict", "kahler_einstein", "dim", "num_rays",
            "degrees", "anticanonical_degree", "mu_TX",
            "max_subsheaf_slope", "witness", "barycenter"]
    recs = [{k: r[k] for k in keep if k in r} for r in rt]
    path = os.path.join(ANC, "root_twist.json.gz")
    write_gzip_json(path, recs)
    print(f"root-twist family: {len(recs)} cases -> root_twist.json.gz "
          f"[{os.path.getsize(path)//1024} KB]")

    with open(os.path.join(OUT, "root_twist_classes.json")) as f:
        cl = json.load(f)
    assert len(cl) == 9 and all(c["verdict"] == "stable" for c in cl)
    assert all(c["kahler_einstein"] ==
               (c["negation_invariant"] or c["rotation_invariant"])
               for c in cl)
    path = os.path.join(ANC, "root_twist_classes.json.gz")
    write_gzip_json(path, cl)
    print(f"zero-sum classes: {len(cl)} classes -> root_twist_classes.json.gz "
          f"[{os.path.getsize(path)//1024} KB]")


README = """Ancillary files for
"Slope stability of tangent bundles of smooth toric Fano varieties"

FILES
  stability_3d.json.gz   all   18 smooth toric Fano 3-folds
  stability_4d.json.gz   all  124 smooth toric Fano 4-folds
  stability_5d.json.gz   all  866 smooth toric Fano 5-folds
  stability_6d.json.gz   all 7622 smooth toric Fano 6-folds
  root_twist.json.gz     the root-twisted dP3-fibrations of Section 7
                         (dimensions 4-8; 7 and 8 beyond the complete census
                         in the paper)
  root_twist_classes.json.gz
                         every zero-sum class of nonzero root twists up
                         to the hexagon's dihedral symmetry (nine
                         classes, dimensions 4-8), with stability,
                         Kahler-Einstein status, and the
                         negation/rotation invariance of the multiset
                         (the table of Section 7)

Each stability file is a gzipped JSON array, one record per variety:

  id                    polyDB identifier (collection
                        Polytopes.Lattice.SmoothReflexive; data of
                        Paffenholz computed by Obro's algorithm)
  picard_rank
  rays                  primitive ray generators of the fan = vertices of
                        the smooth Fano polytope (fan convention; the polyDB
                        source stores the dual polytope -- see Section 3.2)
  degrees               anticanonical degrees deg_{-K} D_rho = normalized
                        volumes of the dual facets, aligned with 'rays'
  anticanonical_degree  (-K)^n  (= sum of 'degrees'; also n! vol of the
                        moment polytope, with both sums formed from the same
                        exact facet triangulations; Section 3.1)
  mu                    slope mu(T_X) = (-K)^n / n as [numerator, denominator]
  max_slope             max over ray-spanned proper subspaces V of
                        mu(F_V) = (sum of degrees of rays in V)/dim V,
                        as [numerator, denominator]
  witness               one subspace attaining max_slope:
                        {dim, rays_in_V (the rays lying in V), deg, slope}.
                        For unstable varieties it is the canonical first
                        Harder-Narasimhan subspace, formed as the sum of all
                        maximal-slope subspaces; for strictly semistable ones
                        it exhibits a slope equality.
  verdict               "stable" | "strictly semistable" | "unstable"
                        (Klyachko criterion, Proposition 2.1 of the paper;
                        all arithmetic exact over Z and Q)
  barycenter            exact barycenter of the moment polytope,
                        list of [numerator, denominator]
  kahler_einstein       true iff barycenter = 0 (Wang-Zhu / Mabuchi)
  polystable            (strictly semistable records only) whether T_X
                        splits equivariantly into stable equal-slope
                        summands, i.e. whether T_X admits a Hermitian-
                        Einstein metric with respect to a Kahler form
                        representing c_1(X) (Proposition 5.2)
  splitting             (polystable records only) a certifying partition of
                        ray indices (0-based into 'rays') whose spans give
                        the equivariant direct-sum decomposition

root_twist.json.gz has the same conventions for the family X(t_1..t_k) of
Definition 7.1, with 'label' encoding the twist tuple.

CONSISTENCY TALLIES (must match the paper)
  dim  total  stable  strictly-ss  unstable  polystable-among-ss  HE-total
   3     18       4         7          7            4                 8
   4    124      26        24         74           14                40
   5    866     109        82        675           52               161
   6   7622     636       366       6620          259               895
  Every Kahler-Einstein variety is stable or polystable (0 unstable KE).
  Barycenters and slopes are different functions of the same exact polytope
  data; comparison with the published low-dimensional KE lists supplies an
  independent end-to-end check.

PROVENANCE
  Generated by src/make_ancillary.py from the exact sweep outputs; the
  generator asserts every tally above against the paper's tables before
  writing. Verdicts for dimensions 3-6 computed by the exact harness
  (src/toric_stability.py, src/sweep.py, src/polystability.py); the input
  smooth-Fano conditions are rechecked by src/validate_input_data.py. The
  dimension-7/8 entries of root_twist.json.gz are independently reproduced
  by the exact rational hexagon formulas in src/prism_check.py.  The hybrid
  accelerated implementation (src/toric_fast.py) agrees with the reference
  implementation wherever both run and supplies an auxiliary cross-check.
  Validation ladder in Section 3 of the paper.
"""


def main():
    os.makedirs(ANC, exist_ok=True)
    for d in (3, 4, 5, 6):
        build_dim(d)
    build_root_twist()
    with open(os.path.join(ANC, "README.txt"), "w") as f:
        f.write(README)
    print("wrote anc/README.txt")
    print("ANCILLARY PACKAGE COMPLETE -- all paper tallies asserted.")


if __name__ == "__main__":
    main()
