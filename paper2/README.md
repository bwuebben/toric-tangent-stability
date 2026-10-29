# Tangent stability on toric Fano blowups and bundles

**Bernd Johannes Wuebben — official active version of Paper 2.**

Original version: September 20, 2026. This version: October 29, 2026.

[Read the paper (PDF, 75 pages)](paper/main.pdf) ·
[LaTeX source](paper/main.tex) · [Exact proof package](paper/anc/README.txt) ·
[Zenodo:23774632](https://zenodo.org/records/23774632) ·
[Abstract and citation](https://bwuebben.github.io/papers/toric-tangent-blowups-bundles/)

This is the second paper in the repository's programme on anticanonical
slope stability of tangent bundles of smooth toric Fano varieties. It develops
uniform geometric results beyond the first paper's classification and
root-twisted families.

## Principal results

- **Blowup trees.** For the specified codimension-two blowups of products of
  projective spaces, every tangent-subspace inequality reduces to one normalized
  divisor-degree inequality per vertex. Exceptional contributions are expressed
  as ordinary mixed intersections on the blowup centres. Every tree of maximum
  degree three is stable when nonleaf factor dimensions equal their valences and
  leaf dimensions are arbitrary positive integers. A saturated degree-four vertex
  surrounded by four saturated degree-four vertices forces instability. Other
  degree-four and degree-five arrangements remain stable.
- **Sharp bundle maximum.** Among the specified locally trivial toric bundles
  with degree-six del Pezzo-product fibres and projective-space-product bases,
  allowing arbitrary integral twists, the maximal stable Picard number is
  `floor(5n/4)+1` for every `n >= 4`. Explicit families attain it.
- **Abundance.** Every fixed subcubic tree with an edge yields stable maximizers
  after sufficiently long independent edge subdivisions. For every integer q,
  all sufficiently large j admit at least q pairwise nonisomorphic toric fans of
  stable extremizers in dimension 4j with Picard number 5j+1.

The two subdivision statements must be distinguished: long subdivisions
stabilize the displayed **bundle** families, while subdivisions can destabilize
members of the **blowup** class. Outside the tree criterion, an inequality for
the specified rank-one foliation does not by itself establish full stability.
The proposed unrestricted Picard bound `rho <= 4n/3+1` is **open**.

## Proofs and computation

The bundle existence proofs combine finite exact rational inequalities with
uniform analytic error bounds and induction to prove entire infinite families.
They are computer assisted. The finite comparisons, matrices and error estimates
are printed in the proof appendices. The blowup-tree arguments use uniform
identities, integral comparisons and induction with explicit polynomial arithmetic.

The [Python programs](paper/anc/code/), [exact inputs and outputs](paper/anc/data/),
and [rational certificates](paper/anc/certificates/) are in `paper/anc/`.
The same package will accompany the arXiv version as ancillary files; Paper 2
is not yet on arXiv.

The theorem-verification commands below require Python 3.11 or newer and only
its standard library. The package also includes the classification screen and
certificate-search sources cited in the paper; these optional searches require
NumPy and SciPy. The [package README](paper/anc/README.txt) gives their commands,
input files and dependencies.
From that directory, quick exact checks include:

```sh
python3 -S code/mixed_star_table.py --check-json data/mixed_star_table/signs.json
python3 -S code/mixed_tree_arithmetic.py --check-json data/mixed_tree_arithmetic/signs.json
python3 -S code/two_junction_gluing.py
python3 -S code/subcubic_thresholds.py
python3 -S code/verify_certificate.py certificates/*.json
```

The longer path, attachment and independent-three-arm computations are also
included; their commands are in [the package README](paper/anc/README.txt).
There are 93 stored rational basis decompositions. Their verifier takes each
file's degree vector as input; the geometric integral formulas provide separate
degree derivations.

The package includes the completed classification slices: all 1,277 sevenfolds
with Picard number at least 9, all 1,582 eightfolds with Picard number at least
11, and the 29 stable boundary examples with Picard number n+2 in dimensions
4–7. An equality witness excludes stability but does not prove semistability.
The new complete Picard-number-10 eightfold census certifies all 6,281 cases:
**64 stable, 692 strictly semistable and 5,525 unstable**. All 692 strictly
semistable fans are nontrivial integral products, so none admits a stable
polarization. The exact exchange graph recovers every tight flat. See
[the census instructions](paper/anc/README-census8.txt); the standalone verifier
needs NumPy and no SciPy. Classification completeness and nonrepetition are
inputs from the cited classification. No dimension-nine census is claimed.

The cycle table generator explicitly cited in the manuscript is
[paper/anc/code/compute06_cycle_tables.py](paper/anc/code/compute06_cycle_tables.py),
with the earlier repository wrapper retained at `figures/code/compute06_cycle_tables.py`.
Both import the portable exact module;
the mathematical computation is unchanged. From `paper/anc/`, run
`python3 -S code/compute06_cycle_tables.py` to write `output/cycle_seeds.tex`.

## Build and version

```sh
cd paper
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

The checked-in PDF and LaTeX sources are the official active version. The
printed original date remains September 20, 2026, with a separate October 29,
2026 version date. This revision adds the complete eightfold census and its
exact geometric, degree, slope and product certificates. Earlier versions
remain in Git history. `SHA256SUMS` records the distributed files; `paper/anc/SHA256SUMS` separately records the
portable mathematical inputs. The historical research draft and unfinished
unrestricted-bound programme are not part of this public manuscript package.

The revised PDF is archived in the existing Zenodo version series. The original
publication date and the current version date are recorded separately.

## Citation

```bibtex
@misc{Wuebben2026ToricBlowupsBundles,
  author = {Wuebben, Bernd Johannes},
  title  = {Tangent stability on toric Fano blowups and bundles},
  year   = {2026},
  note   = {Preprint, original September 20, 2026; this version October 29, 2026},
  doi    = {10.5281/zenodo.23774632},
  url    = {https://bwuebben.github.io/papers/toric-tangent-blowups-bundles/}
}
```

The repository [license](../LICENSE) applies: manuscript and mathematical data
under CC BY 4.0, source code under MIT. Polytope data retain the attribution in
the ancillary README.
