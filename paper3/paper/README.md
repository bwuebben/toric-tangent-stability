# Manuscript and exact supplement

*Local stability chambers for toric tangent bundles*, Bernd Johannes
Wuebben, 26 October 2026, 17 pages.

Read [the PDF](main.pdf) or [the paper guide](../README.md).

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
python3 reproducibility/verify.py
```

The bibliography is included in `main.tex`. The generated inputs
`classification_table.tex` and `ray_data.tex` contain the finite chamber
classification and all 31 ray configurations, degrees, and equality blocks.
All files needed to build are included here.

The [self-contained supplement](reproducibility/README.md) explains the
data conventions, input provenance, and exact scope of the verifier.
The general and all-dimensional results have geometric proofs in the paper;
the finite verifier does not recertify the original exhaustive census.

Rights are stated in [the paper guide](../README.md#rights).
