# Exact chamber computations

Supplement to *Local stability chambers for toric tangent bundles*.

## Run

This directory is self-contained. Copy it anywhere and run:

```bash
python3 verify.py
```

Python 3.10 or later is sufficient; no additional packages are required.
Run without `-O`, since mathematical checks use assertions. A failed
assertion exits nonzero. The computation takes approximately six seconds
on the development machine; timings vary by machine.

The verifier reconstructs all facets by exhaustive supporting equations,
checks unimodularity, computes divisor degrees from exact triangulations
of polar facets, and enumerates flats by intersections of ray-spanned
hyperplanes. It checks all 4,890 proper flats of the 31 indecomposable
cases, their equality blocks, the cone condition for subbundles, numerical
independence of the block divisors, and rational strict curve directions.
It then independently reconstructs the matroid components of all 479
selected configurations by enumerating minimal dependent subsets, and
checks the two projective-bundle fan identifications.

## Files and conventions

- `chambers.json`: 31 complete records, the 448 product partitions, counts
  by dimension, and SHA-256 hashes of eight original source files.
- `selected_rays.json`: ray configurations of all 479 selected cases.
- `projective_bundle_family.json`: finite intersection checks and the
  integral changes of basis identifying the two census examples.
- `verify.py`: standard-library exact verifier. The family data's finite
  intersection checks are supplementary examples; this verifier checks
  its two lattice identifications. The all-dimensional theorem is proved
  in the manuscript.

Ray indices start at zero. Each ray in the JSON is a row vector; these
become columns of the printed matrices. Degrees are ordinary intersection
numbers, without factorial normalization. Fractions in curve directions
are strings parsed exactly by `fractions.Fraction`.

`strict_slope_interval_upper_bound` is an upper bound for strict slope
inequalities along the specified line in curve-class space. It is **not**
a bound on how far that line remains the power of an ample divisor.
The local existence of corresponding polarizations uses the Hodge-index
argument in the manuscript.

## Imported classification

Selection of the 479 strictly semistable cases from the 8,630-variety
census is imported from the preceding paper. This package verifies the
complete local slope data for the 31 indecomposable cases and the matroid
decompositions for all 479 selected cases. The 448 products' semistability
labels and smooth complete fan assumptions remain imported. The package
does not establish census completeness or exclude additional strictly
semistable cases among the other records.

If the original `complex_vb` source tree is available, compare the package
with it by supplying its directory:

```bash
python3 verify.py --source-root /path/to/complex_vb
```

This optional check verifies all eight recorded input hashes, the selected
identifiers, and every copied ray configuration. Even this comparison
does not rerun the original full stability census.

In the source project, `src/build_chamber_data.py` regenerates the 31
records using equivariant localization and closure enumeration. It
requires SciPy for facet proposals, which are subsequently certified by
exact supporting equations and ridge adjacency. The independent verifier
here uses different algorithms and requires no floating-point library.
