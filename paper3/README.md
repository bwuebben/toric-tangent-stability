# Local stability chambers for toric tangent bundles

Bernd Johannes Wuebben · 26 October 2026 · 17 pages.

[Current version (PDF)](paper/main.pdf) · [LaTeX source](paper/main.tex) ·
[Exact computations](paper/reproducibility/README.md) ·
[Zenodo:23681809](https://zenodo.org/records/23681809) ·
[Abstract and citation](https://bwuebben.github.io/papers/toric-local-stability-chambers/)

For the Fano projective bundles
$X_{p,q}=\mathbb P_{(\mathbb P^1)^{p+q}}(\mathcal O\oplus\mathcal O(1^p,-1^q))$,
with nonnegative integers $p,q$ and $p+q\ge1$, the tangent bundle is
anticanonically stable when $p=q$, strictly semistable when $|p-q|=1$,
and unstable otherwise. For $X_{k,k-1}$, the local semistable polarization
locus at the anticanonical class is real analytically equivalent to
$[0,\infty)^k\times\mathbb R^k$, with stable locus
$(0,\infty)^k\times\mathbb R^k$. The number of independent local chamber
boundaries is therefore unbounded.

A general criterion describes local chambers by the order ideals of a finite
poset determined by the equal-slope ray subspaces. Under the stated maximal-cone condition,
nearby stable polarizations exist precisely when the variety is not a
nontrivial toric product. The criterion applies to all 31 indecomposable
strictly semistable cases in the existing dimension-three-through-six
census. Their chambers have one, two, or three independent boundaries.
The paper also determines locally free, pairwise nonisomorphic stable
Jordan–Hölder factors in these cases and in the projective-bundle family.

## Build and reproduce

From `paper3/paper/`:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
python3 reproducibility/verify.py
```

The finite verifier requires Python 3.10 or later and only its standard
library. It reconstructs the local slope data of the 31 indecomposable
examples and the product partitions of all 479 selected configurations.
The exhaustive anticanonical census and the selection of those 479 cases
are imported from [Paper 1](../README.md#paper-1-classification-and-root-twisted-families).
The all-dimensional theorem is proved in the manuscript. See the
[supplement guide](paper/reproducibility/README.md) for the precise scope.

The two generated LaTeX inputs are included, so the manuscript builds
without the research repository. `SHA256SUMS` records the distribution files.
This is a separate preprint; arXiv submission is pending.

## Rights

Copyright © 2026 Bernd Johannes Wuebben. All rights reserved for the
manuscript, code, and original mathematical data in `paper3/`. Existing
attributions and rights in the imported census data are retained. The
licenses for Papers 1 and 2 do not grant additional rights in Paper 3.

## Citation

```bibtex
@misc{Wuebben2026LocalToricChambers,
  author = {Wuebben, Bernd Johannes},
  title = {Local stability chambers for toric tangent bundles},
  year = {2026},
  doi = {10.5281/zenodo.23681809},
  url = {https://doi.org/10.5281/zenodo.23681809}
}
```
