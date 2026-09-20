Ancillary files for
"Tangent stability on toric Fano blowups and bundles".

screen_7d.json.gz: all 1,277 smooth toric Fano sevenfolds with Picard number >=9.
screen_8d.json.gz: all 1,582 smooth toric Fano eightfolds with Picard number >=11.
The eightfolds of Picard number 10 and dimension-nine stability cases are not included.
Sources: Paffenholz's distribution of Obro's classification,
https://polymake.org/polytopes/paffenholz/www/fano.html .
Primitive ray generators are integer vectors. Divisor degrees are exact integers in the
same order. A witness lists a proper ray-spanned subspace of slope >=(-K)^n/n.
An equality witness excludes stability but does not establish semistability.

boundary_stable.json.gz: exactly 29 stable varieties at Picard number n+2,
dimensions 4 through 7, with ray data and structural annotations. The stable eightfold
at Picard number n+3 is excluded from this boundary list.

family.json.gz: the recorded computed twisted del Pezzo bundle cases. This auxiliary
file is not a claim that all relevant families in every dimension have been enumerated.

code/: exact rational programs used for the path, cycle, end-attachment, and
branching-tree theorems. From this directory:
  python3 code/path_all_j.py
  python3 code/path_ends.py
  python3 code/research_graph_targets.py
  python3 code/end_attachment_proof.py --uniform
  python3 code/balanced_spider_proof.py
  python3 code/unequal_spider_proof.py --uniform
  python3 code/two_junction_gluing.py
  python3 code/subcubic_thresholds.py
  python3 -S code/mixed_star_table.py --check-json data/mixed_star_table/signs.json
  python3 -S code/mixed_tree_arithmetic.py --check-json data/mixed_tree_arithmetic/signs.json
The first computes the finite inputs to the equal-end path proof. The second also
runs independent graph-degree, rank, and exhaustive-enumeration checks for all end types.
These programs need Python 3.11 or newer and its standard library only. Their finite checks combine with
the contraction and cycle-deletion arguments printed in the paper to prove all lengths.
The corrected uniform ray-weight bound is 27 and the total error is <6e-17.

The end-attachment proof uses B=243, finite j=2..44 and four comparison lengths;
its common certified margin over the specified invariant configurations is >0.02833. The balanced-spider proof uses B=81,
finite L=2..43 and four comparison lengths; its common margin over the specified invariant configurations is >0.0876.
The unequal-arm theorem extends stability to every triple of arm lengths >=2.
Its proof uses q=16, P=40, 7,140 exact mixed comparisons and a summed error
bound <0.002217. The resulting common margin over configurations invariant
under the local fan symmetries is >0.02698. Arm permutations are not used.
The general-tree theorem builds on this margin by a boundary-rank inequality,
nonnegative limiting transitions and induction on the number of trivalent
vertices. It permits independently varying subdivision lengths above an
explicit threshold depending only on their number; the threshold is O(log N).
The gluing command verifies the limiting weight enclosures, 24 transitions,
12 joining choices and the earlier bridge-33 error bound. The threshold
command reproduces the general theorem's displayed sufficient lengths.
Their outputs are in data/research_two_junction_uniform/constants.json and
data/subcubic_thresholds/thresholds.json. No census of tree shapes is used.
The two-junction corollary with bridge >=27 uses the general error formula
at N=2,q=13, with the inherited three-arm margin and four outer arms >=2.
The portable unequal-arm copy changes only the balanced module's import name.
The portable gluing copy changes only the balanced and unequal module names;
the subdivision threshold program is copied unchanged.
The E1 portable copy differs from the reviewed research module only in its optional
matroid_certificate import and the no-argument entry point. It imports the identical
exact rref/in_span/rank functions directly from toric_stability, avoiding NumPy/SciPy.
Optional numerical basis search remains in the original research source and is not
required to verify the stored certificates or any all-length computation.

The mixed-star command verifies all fourteen fractions in the star-transition
table directly from the printed finite volume formula and its derivative.
It is standalone, is copied unchanged, and uses only exact rational arithmetic.
Its check ends at m=14; the proof for all m>=15 is the analytic estimate in
the manuscript. The tree criterion, all-length mixed-blowup paths and
saturated-vertex comparison do not depend on a finite graph computation.

The mixed-tree arithmetic command reproduces the printed reference polynomials,
closure inequalities and root degree for all trees of maximum degree three;
the degree-four and degree-five local obstructions; and the endpoint signs,
eighteen linear-functional values and six adjacent pairs in the two-centre
proof. It is standalone and copied unchanged. Arbitrary tree sizes follow from
the ratio-order induction, arbitrary outside branches from the comparison
theorems, and arbitrary connecting lengths from positivity of polynomial
coefficients. Finite graph diagnostics are not premises of these conclusions.

data/: frozen exact comparison outputs. The graph and new-family commands write output/;
compare the resulting JSON objects exactly with data/. Elapsed time appears only
in stdout, not in the frozen JSON. The graph copy changes only its stored
spider222 certificate path to the package layout. These finite
computations supplement the printed geometric, rank, contraction and deletion proofs.

certificates/: 93 stored positive rational combinations of incidence vectors of ray bases,
including all 20 short end-attachment examples.
Check every basis, coefficient, identity and exchange digraph with:
  python3 code/verify_certificate.py certificates/*.json
The verifier takes the normalized degrees in each file as input; it does not independently
recompute them. The transfer formulas in the paper and code compute path degrees.
P8D_389596.json uses degrees obtained directly from the classification rays.

SHA256SUMS records all package contents other than itself.
