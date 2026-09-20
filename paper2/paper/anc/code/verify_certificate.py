"""Standalone exact verifier for output/certificates/*.json (no imports from this project).
Checks, in rational arithmetic: delta > 0 and sum delta = n = rank of the rays; every B_k is a
basis; lambda_k > 0, sum lambda_k = 1, sum_k lambda_k 1_{B_k} = delta; the exchange digraph
(e -> e' iff for some k: e not in B_k, e' in B_k, B_k - e' + e a basis) is strongly connected.
By Lemmas A and B of matroid_certificate.py this proves stability, given that delta is the vector
of normalized degrees n d_v / D (its source is recorded in the file; it is not re-derived here).
Usage: verify_certificate.py file.json [...]"""
import sys, json
from fractions import Fraction as Fr

def solve_in_basis(rows, targets):
    """rows: n independent vectors.  Returns, for each target t, the coefficients c with t = sum c_k rows[k];
    raises if the rows are dependent."""
    n = len(rows); M = [[Fr(rows[k][i]) for k in range(n)] + [Fr(t[i]) for t in targets] for i in range(n)]
    for c in range(n):
        p = next((i for i in range(c, n) if M[i][c] != 0), None)
        if p is None: raise ValueError("not a basis")
        M[c], M[p] = M[p], M[c]; inv = 1 / M[c][c]; M[c] = [x*inv for x in M[c]]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c]; M[i] = [a - f*b for a, b in zip(M[i], M[c])]
    return [[M[k][n + t] for k in range(n)] for t in range(len(targets))]

def verify(path):
    C = json.load(open(path)); rays = [tuple(v) for v in C['rays']]; E, n = len(rays), len(rays[0])
    delta = [Fr(x) for x in C['delta']]; lam = [Fr(x) for x in C['lambdas']]; bases = [tuple(B) for B in C['bases']]
    assert all(x > 0 for x in delta) and sum(delta) == n, "delta"
    assert all(x > 0 for x in lam) and sum(lam) == 1, "lambda"
    assert all(sum(l for B, l in zip(bases, lam) if e in B) == delta[e] for e in range(E)), "decomposition"
    arcs = [set() for _ in range(E)]
    for B in bases:
        assert len(set(B)) == n, "basis size"
        outside = [e for e in range(E) if e not in B]
        coeffs = solve_in_basis([rays[q] for q in B], [rays[e] for e in outside])       # raises if B is dependent
        for e, cf in zip(outside, coeffs):
            arcs[e].update(B[k] for k in range(n) if cf[k] != 0)
    def reach(adj):
        seen = {0}; stack = [0]
        while stack:
            x = stack.pop()
            for y in adj[x]:
                if y not in seen: seen.add(y); stack.append(y)
        return len(seen) == E
    rev = [set() for _ in range(E)]
    for x in range(E):
        for y in arcs[x]: rev[y].add(x)
    assert reach(arcs) and reach(rev), "exchange digraph not strongly connected"
    return C['name'], n, E, len(bases)

if __name__ == '__main__':
    for p in sys.argv[1:]:
        name, n, E, nb = verify(p)
        print(f"VERIFIED  {name}: n = {n}, {E} rays, {nb} bases with positive rational weights summing to delta; exchange digraph strongly connected -> stable")
