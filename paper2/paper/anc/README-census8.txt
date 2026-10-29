Exact certificates for the Picard-number-ten eightfold stratum
29 October 2026

Contents
  code/verify_census8.py
  data/census8_certificates.json.gz
  data/smooth_toric_fano_8d_rho_ge_10.json (the existing classification input)
  requirements-census8.txt
  SHA256SUMS-census8

Run from this ancillary directory with Python 3.11 or newer:
  python3 -m pip install -r requirements-census8.txt
  OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -B code/verify_census8.py --workers 2 --report output/census8_verification.json

Check the distributed file hashes with:
  shasum -a 256 -c SHA256SUMS-census8

The paths are resolved relative to the verifier, so the command also works
from another working directory. --input and --certificates accept alternate
locations for the two data files. The input must retain the supplied checksum.
The verifier reads no other modules, results, reports, or classification code.
SciPy is not required. --pilot checks twelve deterministic positions and does
not certify the full population. The complete command has no cache and checks
every record on every run. A failed check produces a nonzero exit status.

Scope and expected output
The raw file is Paffenholz's distribution of Obro's smooth Fano classification,
restricted to dimension eight and Picard number at least ten:
  https://polymake.org/polytopes/paffenholz/www/fano.html
The verifier selects exactly the 6,281 records with eighteen rays, hence Picard
number ten. Their completeness and pairwise nonisomorphism are imported from
that classification; the database is not independently regenerated here.
The higher-Picard-number records in the same raw file are outside this check.

All 6,281 records are certified: 64 stable, 692 strictly semistable, and 5,525
unstable. Verification covers 2,443,089 maximal cones, 113,058 divisor degrees,
and 12,329 linear bases. The 692 strictly semistable records are also checked
as nontrivial toric products, including their full Cartesian-product fans,
integral lattice decompositions, complete equality lists, and independence
of the equality-block divisor classes modulo principal divisors. There are
4,208 proper equalities. Every maximal cone contains exactly the rank-many
rays of each proper equality flat: no failure of this condition occurs.
The number of finest product components has distribution
  2:464, 3:170, 4:44, 5:10, 6:3, 7:1.
These conclusions apply only to this supplied stratum.

Verification method
Each proposed maximal cone is an eight-ray unimodular basis with a strict
integral supporting equation. Ridges occur twice, facet adjacency is connected,
and every ray occurs, certifying the complete fan. NumPy inverses are proposals
only: multiplication by the original matrix over Python integers certifies
each accepted inverse. Exact elimination supplies a fallback.

Rational localization recomputes every divisor degree using a deterministic
integral vector different from the stored vector. The verifier also checks
positivity, integrality, the anticanonical sum, and principal-divisor relations.
An unstable record has a closed proper ray flat with positive normalized slope
excess. A semistable record has a strictly positive rational convex combination
of ray bases reconstructing the normalized degrees. Exact basis coordinates
reconstruct its exchange graph. Strong connectedness distinguishes stability;
a strictly semistable record additionally has a tight sink-component witness.
An equality witness by itself is never used to prove semistability.

Product components are recomputed from the full rational nullspace. A maximal
cone supplies a unimodular basis separating their integral lattices. Projection
cardinalities certify that the fan is the full Cartesian product. All proper
outgoing-closed unions of exchange strong components are enumerated and checked
as equality flats. Rational elimination computes equality-block divisor classes
modulo principal divisors. block_ranks in the bundle are positive quotient
dimensions (sums of normalized degrees), which need not equal the dimensions
of the raw spans of the individual blocks. The report prints both quantities.

The compressed bundle contains logical certificate data and no run timings.
Its JSON has records sorted by identifier. Each strictly semistable record
includes a structure object. The detailed verification report is deterministic
and records the raw-input, bundle, verifier, and canonical-record SHA-256 hashes.
The SHA-256 of the unchanged raw input is
  fec018c5d465ebe8ec664b1a698710df7be6cdfef5967c50d1646e57237e58ee
