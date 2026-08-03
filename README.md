# Slope stability of tangent bundles of smooth toric Fano varieties

Reproducibility and verification companion to the paper

> **Slope stability of tangent bundles of smooth toric Fano varieties**
> Bernd Johannes Wuebben, 2026.
> arXiv: *(to be added on posting)*

The paper determines the slope stability, with respect to the anticanonical
polarization, of the tangent bundle of **every smooth toric Fano variety of
dimension at most six** — 18, 124, 866 and 7622 varieties in dimensions 3, 4,
5 and 6 — by exact evaluation of Klyachko's criterion. It further
decides polystability for every strictly semistable case (the
Hermitian–Einstein census): the connected components of the ray matroid give
the finest equivariant splitting of the tangent bundle, so a strictly
semistable smooth toric Fano has polystable tangent bundle exactly when it is
a nontrivial product of factors with stable tangent bundles. It compares the
census with the Wang–Zhu Kähler–Einstein criterion and introduces the
*root-twist* construction: toric dP₃-fibrations over products of projective
lines twisted by roots of A₂. Every nonempty multiset of nonzero root twists
with vanishing sum yields a stable tangent bundle, in every dimension, and
among these the variety is Kähler–Einstein exactly when the multiset is
invariant under negation or the order-three rotation of the root hexagon;
a separate unbalanced family is stable but not Kähler–Einstein in every even
dimension. The proofs reduce the slope inequalities to two-dimensional
hexagon integrals.

Every verdict in the paper is reproducible from the ray data in this
repository, with exact (integer/rational) arithmetic and no external data
dependencies.

## Layout

| Path | Contents |
|---|---|
| `paper/` | The paper: `main.tex` (self-contained, bibliography included), `toric-tangent-stability.pdf`, and `anc/` — the ancillary data package (per-variety degrees, slopes, maximal-slope subspace and verdict for all 8630 varieties, plus the root-twist family and the nine zero-sum twist classes, gzipped, with its own `README.txt`). |
| `src/` | The computation and verification code (Python, standard library + NumPy/SciPy; see below). |
| `data/` | Vertex data of the smooth Fano polytopes (primitive ray generators of the fans) in dimensions 3–6. Each file records its `source` in a `source` field. |
| `results/` | Precomputed verdict tables: `sweep_<n>d.json` (per-variety degrees, slopes, witness subsheaf, verdict, barycenter) and `polystability_<n>d.json` (the polystability refinement), plus `root_twist.json` (the root-twist family, dimensions 4–8) and `root_twist_classes.json` (all nine zero-sum twist classes up to the hexagon's dihedral symmetry, with stability, Kähler–Einstein status, and multiset symmetries). Fully regenerable from `src/` + `data/`. |
| `docs/` | Supplementary mathematics: the Section-7 prism-reduction derivations and the dimension-4 (Picard rank ≤ 3) validation write-up. |

## The method, in one paragraph

Klyachko's classification of equivariant vector bundles makes an equivariant
sheaf on a smooth toric variety into filtration data on the rays; uniqueness of
the Harder–Narasimhan filtration makes the maximal destabilizing subsheaf
equivariant; for the tangent bundle the filtrations are two-step, so
(semi)stability with respect to `-K` reduces to finitely many exact integer
comparisons of slopes over subspaces spanned by subsets of rays. The smooth
toric Fano varieties are completely classified, so the strategy is: a certified
exact sweep over the classification, then pattern, then theorem.

## Reproducing the results

Requirements: Python 3.9+. The core checker `src/toric_stability.py` is
**dependency-free** (pure `fractions`); the accelerated variant
`src/toric_fast.py` and some check scripts additionally need NumPy and SciPy.

```bash
python3 -m venv venv
venv/bin/pip install -r requirements.txt      # numpy, scipy
```

Self-test the exact checker (hand-verifiable examples in dimensions 2–4):

```bash
python3 src/toric_stability.py
```

Re-run a full-dimension sweep (writes a `sweep_<n>d.json` matching `results/`):

```bash
python3 src/sweep.py 3        # dimensions 3, 4, 5; dimension 6 runs sharded
```

Verify the Section-7 claims (exact):

```bash
venv/bin/python src/prism_check.py         # 98/98 cases
venv/bin/python src/inequality_check.py    # claims 1–7
venv/bin/python src/roottwist_classes.py   # all nine zero-sum twist classes
```

Revalidate the input data (reconstructs every fan-polytope facet exactly and
checks the smooth-Fano conditions for all 8,630 inputs):

```bash
python3 src/validate_input_data.py
```

Regenerate the ancillary data package (asserts every tally against the
paper before writing):

```bash
python3 src/make_ancillary.py
```

## Data provenance

The polytope databases in `data/` are the primitive ray generators of the fans
of the smooth toric Fano varieties, obtained from the standard classification
(Watanabe–Watanabe and Batyrev in dimension 3, Batyrev and Sato in dimension 4,
and Øbro's algorithm in dimensions 5–6), as collected in
[**polyDB**](https://polydb.org) (data collection
`Polytopes.Lattice.SmoothReflexive`, Paffenholz). Each JSON file's `source`
field records its origin. Please credit polyDB when reusing the polytope data.

## Building the paper

```bash
cd paper
latexmk -pdf main.tex
```

## License

- **Code** (`src/`): MIT License — see `LICENSE`.
- **Paper** (`paper/`) and **data/results** (`data/`, `results/`, `paper/anc/`):
  Creative Commons Attribution 4.0 International (CC BY 4.0) — see `LICENSE`.

## Citation

```bibtex
@misc{Wuebben2026ToricTangent,
  author = {Wuebben, Bernd Johannes},
  title  = {Slope stability of tangent bundles of smooth toric Fano varieties},
  year   = {2026},
  note   = {arXiv preprint}
}
```
