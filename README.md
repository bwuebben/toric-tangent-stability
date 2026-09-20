# Tangent-bundle stability on smooth toric Fano varieties

Manuscripts, exact computations and data for two related papers.

| Paper | Manuscript and materials |
|---|---|
| **1. Slope stability of tangent bundles of smooth toric Fano varieties** | [Current PDF](paper/toric-tangent-stability.pdf) · [arXiv:2608.20411](https://arxiv.org/abs/2608.20411) · code and data below |
| **2. Tangent stability on toric Fano blowups and bundles** | [PDF](paper2/paper/main.pdf) · [Guide, source and exact proof package](paper2/README.md) · preprint, 20 September 2026 |

## Paper 1: classification and root-twisted families

> **Slope stability of tangent bundles of smooth toric Fano varieties**
> Bernd Johannes Wuebben, 2026.
> [arXiv:2608.20411](https://arxiv.org/abs/2608.20411) (math.DG, cross-listed math.AG)

**Current manuscript:** [September 20, 2026 revision](paper/toric-tangent-stability.pdf)
(original version: August 13, 2026). The arXiv replacement is pending.

The paper constructs smooth toric Fano n-folds of Picard number n + 2 whose
tangent bundles are slope-stable with respect to the anticanonical polarization,
for every n ≥ 4. These are toric fibrations with fibre the del Pezzo surface of
degree six over products of projective lines, parametrized by multisets of roots
of A₂. Every nonempty multiset of nonzero roots with vanishing sum yields a stable
tangent bundle. Among these multisets, the variety is Kähler–Einstein exactly
when the multiset is invariant under negation or the order-three rotation of the
root hexagon. A separate family with nonzero twist sum is stable but not
Kähler–Einstein in every even dimension at least four. The proofs reduce slope
inequalities and barycenter computations to integrals over the moment hexagon.

For a smooth toric Fano variety with strictly semistable tangent bundle,
polystability is equivalent to decomposition as a nontrivial product of smooth
toric Fano varieties with anticanonically stable tangent bundles. Exact evaluation
of Klyachko's criterion extends the stability classification to dimensions five
and six, recovering the classifications of Steffens and Reynolds in dimensions
three and four. Together with the product criterion, this determines polystability
for all **8,630 varieties in dimensions three through six**, and hence the existence
of Hermitian–Einstein metrics on their tangent bundles with respect to an
anticanonical Kähler form.

Every verdict in the paper is reproducible from the ray data in this
repository, with exact (integer/rational) arithmetic and no external data
dependencies.

## Paper 2: stability mechanisms and abundant extremizers

The second paper proves an all-subspace stability criterion for blowup trees,
uniform stability for trees of maximum degree three with the stated factor
dimensions, and a specified degree-four instability obstruction. An exceptional
mixed-intersection formula identifies the degree of the corresponding toric
foliation; stable degree-four/five examples show why valence alone is insufficient.

For the specified toric bundles with degree-six del Pezzo-product fibres and
projective-space-product bases, it proves the sharp stable Picard maximum
`floor(5n/4)+1` for every `n >= 4`, allowing arbitrary integral twists.
Every fixed subcubic tree with an edge gives further maximizers after sufficiently
long independent subdivisions. The resulting number of torically distinct
extremizers is unbounded. The unrestricted `4n/3+1` upper bound remains open.

The bundle existence proofs use exact rational calculations together with uniform
error estimates and induction. Complete proofs, reproducible code and data are
included in [paper2/](paper2/README.md). This is a separate preprint; it is not an
arXiv replacement for Paper 1.

## Layout

| Path | Contents |
|---|---|
| `paper/` | Paper 1: `main.tex` (self-contained, bibliography included), `toric-tangent-stability.pdf`, and `anc/` — the ancillary data package (per-variety degrees, slopes, maximal-slope subspace and verdict for all 8630 varieties, plus the root-twist family and the nine zero-sum twist classes, gzipped, with its own `README.txt`). |
| `paper2/` | Paper 2 guide, complete manuscript and source, and a portable exact proof package with its own README and checksums. |
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

## Reproducing the Paper 1 results

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

## Building the manuscripts

```bash
cd paper
latexmk -pdf -jobname=toric-tangent-stability main.tex
```

For Paper 2, run `latexmk -pdf main.tex` from `paper2/paper/`. Its portable
proof commands and Python requirements are in [paper2/README.md](paper2/README.md).

## License

- **Code** (`src/`, `paper2/paper/anc/code/`, `paper2/figures/code/`): MIT License — see `LICENSE`.
- **Manuscripts and mathematical data** (`paper/`, `data/`, `results/`, and
  the non-code contents of `paper2/`):
  Creative Commons Attribution 4.0 International (CC BY 4.0) — see `LICENSE`.

## Citation

```bibtex
@misc{Wuebben2026ToricTangent,
  author        = {Wuebben, Bernd Johannes},
  title         = {Slope stability of tangent bundles of smooth toric Fano varieties},
  year          = {2026},
  eprint        = {2608.20411},
  archivePrefix = {arXiv},
  primaryClass  = {math.DG},
  note          = {arXiv:2608.20411}
}
```

For Paper 2, use the separate [citation entry](paper2/README.md#citation).
