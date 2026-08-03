# Validation: tangent-bundle (semi)stability of smooth toric Fano 4-folds, ρ ≤ 3, vs. the literature

**Date:** 2026-07-10.
**Our data:** `output/dim4_rho_le3_for_validation.json` (38 varieties: 1 with ρ=1, 9 with ρ=2,
28 with ρ=3), verdicts from our Klyachko facet-volume computation of slope-(semi)stability of
T_X w.r.t. −K, plus a barycenter-criterion KE flag.

**Result: 38/38 matched, 0 disagreements** — on stability verdicts (vs. Dasgupta–Dey–Khan and
Biswas–Dey–Genc–Poddar) and on the KE flag (vs. Nitta–Saito–Yotsutani / Nakagawa).

## Sources

1. **[DDK]** J. Dasgupta, A. Dey, B. Khan, *Stability of equivariant vector bundles over toric
   varieties*, arXiv:1910.13964, Documenta Math. 26 (2021), DOI 10.4171/dm/785.
   - Verdicts extracted from **Table 1 ("Stability of tangent bundle of toric Fano 4-folds"),
     end of Section 5** (pp. 33–34 of the arXiv PDF). The table covers all 38 varieties with
     ρ ≤ 3 in Batyrev's notation (P⁴; B₁–B₅, C₁–C₄; D₁–D₁₉, E₁–E₃, G₁–G₆), citing
     Proposition 4.1.1, Remark 4.2.6, Corollary 4.2.7 (ρ ≤ 2) and Propositions 5.1.1–5.7.4 (ρ = 3).
   - **Fan data used for matching** comes from DDK's own explicit ray lists in the proofs:
     §5.1 (D₇, D₁₇), §5.2 (D₁, D₆, D₁₈, D₁₉), §5.3 (D₂, D₃, D₅, D₈, D₉, D₁₂, D₁₆), §5.4 (D₁₁),
     §5.5 (D₄, D₁₀; also the 3-fold factor of D₁₄), §5.6 (E₁–E₃ and the B-type fans),
     §5.7.1–5.7.4 (G₁–G₆ and the C-type fans).
2. **[BDGP]** I. Biswas, A. Dey, O. Genc, M. Poddar, *On stability of tangent bundle of toric
   varieties*, arXiv:1808.08701 (Proc. Indian Acad. Sci.). **Theorem 9.3** (Section 9) — second
   witness for ρ ≤ 2: (1) T_{P⁴} stable; (2) T_{B₄}, T_{C₄} polystable; (3) T_{B₅} strictly
   semistable (Theorem 9.2); (4) T_{B₁}, T_{B₂}, T_{B₃} unstable; (5) T_{C₁}, T_{C₂}, T_{C₃}
   unstable. Fully consistent with DDK (see caveat on "polystable" below).
3. **[NSY]** Y. Nitta, S. Saito, N. Yotsutani, *Relative Ding and K-stability of toric Fano
   manifolds in low dimensions*, arXiv:1712.01131 (v5). **Table 3** lists all 124 toric Fano
   4-folds in Batyrev/Sato notation with the Mabuchi constant M_X. For toric Fanos, M_X = 0 ⟺
   vanishing Futaki character ⟺ barycenter of the moment polytope is 0 ⟺ KE (Wang–Zhu). Among
   ρ ≤ 3 (Table 3, rows 1–38), M_X = 0 exactly for **P⁴, B₄, C₄, D₁₃** — this is the KE
   cross-check set (agrees with Nakagawa's KE classification, which BDGP §9 also cites).
4. Batyrev, *On the classification of toric Fano 4-folds*, arXiv:math/9801107 — the naming
   scheme (Section 4); not needed for fan data since DDK reproduces the fans explicitly.

## Matching methodology

For each Batyrev name we built the fan's ray-generator set from DDK's explicit constructions
(and standard product/projective-bundle fans for P⁴, B₁–B₅, C₁–C₄, D₁₃, D₁₄, D₁₅, using the
same conventions DDK uses in §5.5–5.7). Each candidate vertex set was then matched against the
38 polyDB vertex sets by **exact GL(4,ℤ) lattice equivalence** (exhaustive search over ordered
4-tuples of target vertices, exact rational arithmetic; script:
`match.py` in the session scratchpad). The outcome is a **perfect bijection**: every polyDB id
matched exactly one Batyrev name and every name was used exactly once. No matching was decided
by invariants alone; degrees were only used as an a-posteriori consistency check.

**Independent degree check:** DDK reports µ(T_X) = (−K)⁴/4 in each proposition. Every reported
value matches the polyDB `anticanonical_degree` of the matched entry, e.g. D₁₇: µ=101.25 → 405
(F.4D.0085); D₁₉: µ=100 → 400 (F.4D.0113); D₁₈: µ=100 → 400 (F.4D.0019); E₃: µ=107.75 → 431
(F.4D.0103); G₄: µ=104.25 → 417 (F.4D.0041); G₅: µ=101.5 → 406 (F.4D.0042); G₆: µ=100.25 → 401
(F.4D.0108); E₁: 151.25 → 605; G₁: 132.25 → 529; etc. This disambiguates in particular the
degree-degenerate pairs (400: D₁₈/D₁₉; 432: D₁₃/D₁₄/D₁₅/D₁₆; 464: D₉/D₁₀; 496: D₅/D₆;
512: B₄/B₅; 513: C₂/C₃; 560: D₃/D₄), which the lattice matching resolves.

## Comparison table (all 38 entries)

Verdicts: S = stable, SSS = strictly semistable (semistable, not stable), U = unstable.
DDK verdicts from Table 1; KE(lit) from NSY Table 3 (M_X = 0).

| polyDB id | ρ | (−K)⁴ | Batyrev name (DDK description) | DDK verdict | our verdict | KE(lit) | our KE | MATCH |
|---|---|---|---|---|---|---|---|---|
| F.4D.0123 | 1 | 625 | P⁴ | S | S | yes | yes | OK |
| F.4D.0001 | 2 | 800 | B₁ = P(O⊕O(3)) over P³ | U | U | no | no | OK |
| F.4D.0115 | 2 | 640 | B₂ = P(O⊕O(2)) over P³ | U | U | no | no | OK |
| F.4D.0120 | 2 | 544 | B₃ = P(O⊕O(1)) over P³ | U | U | no | no | OK |
| F.4D.0121 | 2 | 512 | B₄ = P¹×P³ | SSS | SSS | yes | yes | OK |
| F.4D.0114 | 2 | 512 | B₅ = P(O³⊕O(1)) over P¹ | SSS | SSS | no | no | OK |
| F.4D.0020 | 2 | 594 | C₁ = P(O²⊕O(2)) over P² | U | U | no | no | OK |
| F.4D.0117 | 2 | 513 | C₂ = P(O²⊕O(1)) over P² | U | U | no | no | OK |
| F.4D.0046 | 2 | 513 | C₃ = P(O⊕O(1)²) over P² | U | U | no | no | OK |
| F.4D.0122 | 2 | 486 | C₄ = P²×P² | SSS | SSS | yes | yes | OK |
| F.4D.0008 | 3 | 592 | D₁ = P(O⊕O(1,2)) over P¹×P² | U | U | no | no | OK |
| F.4D.0009 | 3 | 576 | D₂ = P(O⊕O(0,1)) over B₁ | U | U | no | no | OK |
| F.4D.0032 | 3 | 560 | D₃ = P(O⊕O(1,1)) over B₂ | U | U | no | no | OK |
| F.4D.0013 | 3 | 560 | D₄ = P(O⊕O(0,2)) over B₃ | U | U | no | no | OK |
| F.4D.0018 | 3 | 496 | D₅ = P¹×P(O⊕O(2)) over P² | U | U | no | no | OK |
| F.4D.0105 | 3 | 496 | D₆ = P(O⊕O(1,1)) over P¹×P² | U | U | no | no | OK |
| F.4D.0073 | 3 | 486 | D₇ = P(O⊕O⊕O(1,1)) over P¹×P¹ | U | U | no | no | OK |
| F.4D.0107 | 3 | 480 | D₈ = P(O⊕O(0,1)) over B₂ | U | U | no | no | OK |
| F.4D.0043 | 3 | 464 | D₉ = P(O⊕O(1,0)) over B₂ | U | U | no | no | OK |
| F.4D.0106 | 3 | 464 | D₁₀ = P(O⊕O(0,1)) over B₃ | U | U | no | no | OK |
| F.4D.0093 | 3 | 459 | D₁₁ = P(O⊕O⊕O(0,1)) over H₁ | U | U | no | no | OK |
| F.4D.0116 | 3 | 448 | D₁₂ = P¹×P(O⊕O(1)) over P² | U | U | no | no | OK |
| F.4D.0119 | 3 | 432 | D₁₃ = P¹×P¹×P² | SSS | SSS | yes | yes | OK |
| F.4D.0111 | 3 | 432 | D₁₄ = P¹×P(O²⊕O(1)) over P¹ | SSS | SSS | no | no | OK |
| F.4D.0112 | 3 | 432 | D₁₅ = H₁×P² | SSS | SSS | no | no | OK |
| F.4D.0045 | 3 | 432 | D₁₆ = P(O⊕O(−1,1)) over B₂ | U | U | no | no | OK |
| F.4D.0085 | 3 | 405 | D₁₇ = P(O⊕O(1,0)⊕O(0,1)) over P¹×P¹ | S | S | no | no | OK |
| F.4D.0019 | 3 | 400 | D₁₈ = P(O⊕O(−1,2)) over P¹×P² | U | U | no | no | OK |
| F.4D.0113 | 3 | 400 | D₁₉ = P(O⊕O(−1,1)) over P¹×P² | S | S | no | no | OK |
| F.4D.0000 | 3 | 605 | E₁ = Bl_{P²}(B₂) | U | U | no | no | OK |
| F.4D.0104 | 3 | 489 | E₂ = Bl_{P²}(B₃) | U | U | no | no | OK |
| F.4D.0103 | 3 | 431 | E₃ = Bl_{P²}(B₄) | S | S | no | no | OK |
| F.4D.0003 | 3 | 529 | G₁ | U | U | no | no | OK |
| F.4D.0017 | 3 | 450 | G₂ = Bl_{P¹×P¹}(C₂) | U | U | no | no | OK |
| F.4D.0026 | 3 | 433 | G₃ = Bl_{P¹}(C₃) | U | U | no | no | OK |
| F.4D.0041 | 3 | 417 | G₄ = Bl_{H₁}(C₂) | S | S | no | no | OK |
| F.4D.0042 | 3 | 406 | G₅ = Bl_{P¹×P¹}(C₃) | S | S | no | no | OK |
| F.4D.0108 | 3 | 401 | G₆ = Bl_{P¹×P¹}(C₄) | S | S | no | no | OK |

Here H₁ = Hirzebruch surface F₁, B₁ = P(O⊕O(2)) over P³ etc. as in Batyrev; DDK's Bl_Z(X)
notation names the blow-up centre by its isomorphism type.

## Aggregate consistency

- **ρ = 3 (28 varieties):** DDK (Introduction and Table 1): 6 stable (D₁₇, D₁₉, E₃, G₄, G₅, G₆),
  3 strictly semistable (D₁₃, D₁₄, D₁₅), 19 unstable. Ours: identical, variety by variety.
- **ρ = 2 (9 varieties):** DDK/BDGP: 0 stable, 3 strictly semistable (B₄, B₅, C₄), 6 unstable.
  Ours: identical.
- **ρ = 1:** P⁴ stable (both).
- **KE (ρ ≤ 3):** literature list {P⁴, B₄, C₄, D₁₃} = our KE=true set
  {F.4D.0123, F.4D.0121, F.4D.0122, F.4D.0119}. Exact agreement. (BDGP §9 explicitly notes B₅
  is strictly semistable but *not* KE — our data agrees: F.4D.0114 has KE=false.)

## Disagreements

**None.** All 38 stability verdicts and all 38 KE flags agree with the literature.

## Caveats on notation / terminology mapping

- **"Polystable" (BDGP Thm 9.3(2)) vs. "strictly semistable" (DDK Table 1) for B₄, C₄.** For
  these product varieties T_X is a direct sum of equal-slope stable pieces, hence polystable and
  in particular semistable-but-not-stable. Both statements are consistent with our
  "strictly semistable" verdict; there is no conflict.
- **Sign conventions in P(E).** DDK builds fans for P(O⊕O(D)) etc. with the twist recorded on
  the base ray; since P(E) ≅ P(E⊗L), several (α,β) choices give isomorphic varieties. We used
  DDK's exact printed ray lists wherever given (all D, E, G types), so name assignment follows
  DDK's own identification with Batyrev's [math/9801107, §4] names (which DDK states explicitly
  in each proposition, e.g. "X = D₁, D₆, D₁₈, D₁₉ when (α,β) = (1,2), (1,1), (−1,2), (−1,1)").
- **Degenerate invariants.** Several (ρ, #rays, degree) triples are shared by 2–4 varieties
  (degrees 400, 432, 464, 496, 512, 513, 560). All such pairs were disambiguated by exact
  GL(4,ℤ) equivalence of vertex sets, not by invariants.
- **G₃'s blow-up centre.** DDK's Table 1 prints G₃ = Bl_{P¹}(C₃) while the proof (Prop. 5.7.3)
  blows up V(τ) for the 3-dimensional-quotient τ = Cone(v₁,v₂,e₀) ⊂ Δ′(C₃), i.e. a curve; the
  printed fan (rays incl. u_τ = (1,1,−1,−1)) is what we matched against, so the name assignment
  is unaffected by how the centre is labelled.
- **NSY notation order.** NSY Table 3 rows 1–38 enumerate exactly P⁴, B₁–B₅, C₁–C₄, E₁–E₃,
  D₁–D₁₉, G₁–G₆ (Batyrev names), confirming that these 38 names exhaust ρ ≤ 3 — the same count
  as our JSON.
- NSY's Table 3 concerns *relative Ding stability* (Mabuchi constants), **not** slope stability
  of T_X; we used it only for the KE flag via M_X = 0 (vanishing Futaki ⟺ KE for toric Fanos,
  Wang–Zhu) and for the notation census. E.g. G₄–G₆ are relatively Ding *unstable* yet have
  slope-*stable* tangent bundle — different notions, no contradiction.

## Bottom line

Our Klyachko-criterion computation reproduces the published classification of tangent-bundle
(semi)stability for **all 38** smooth toric Fano 4-folds of Picard rank ≤ 3, and our barycenter
KE flags reproduce the known KE set in this range. The validation set passes completely; the
ρ ≥ 4 verdicts (86 varieties) are new territory not covered by DDK/BDGP.
