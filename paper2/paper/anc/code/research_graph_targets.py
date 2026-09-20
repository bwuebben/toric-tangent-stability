"""Bounded exact checks for the September 2026 graph research programme.

No pre-existing manuscript or proof modules are modified.  The cycle argument
and the graph rank formula are written in research_notes/.../graph_targets.md.
The default checks use only the standard library and finish in a short run.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import math
from functools import lru_cache
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import path_transfer as P
import tree_family as G

OUT = Path(__file__).resolve().parents[1] / 'output' / 'research_graph_targets'


def trace(A):
    return sum(A[i][i] for i in range(3))


def cycle(j):
    assert j >= 2
    return G.HexGraph([2]*j, [((i+1) % j, i) for i in range(j)], f'cycle{j}')


def cycle_weights(j):
    A = [[F(i == k) for k in range(3)] for i in range(3)]
    for _ in range(j-1):
        A = P.mm(A, P.T)
    Z = trace(P.mm(A, P.T))
    H = {kind: trace(P.mm(P.mm(W, P.KAP), A)) / Z
         for kind, W in P.WH.items()}
    c = trace(P.mm(P.mm(P.W, P.KAP1), A)) / Z
    assert H['r'] == H['-r'] and H['neg'] == H['pos']
    a, b = H['r'], H['neg']
    assert 2*a + 4*b + 3*c == 4
    return Z, a, b, c


def graph_rank(graph, H, C, T):
    """Proposed general graph formula, independently compared with rational rref."""
    adj = [[] for _ in graph.b]
    for i, (p, m) in enumerate(graph.edges):
        if H[i] == 0:
            adj[p].append(m)
            adj[m].append(p)
    remaining = set(range(graph.k))
    complete = 0
    while remaining:
        start = next(iter(remaining))
        stack, comp = [start], {start}
        remaining.remove(start)
        while stack:
            for v in adj[stack.pop()]:
                if v in remaining:
                    remaining.remove(v)
                    comp.add(v)
                    stack.append(v)
        complete += all(C[v] and T[v] for v in comp)
    return sum(H) + sum(b*x for b, x in zip(graph.b, C)) + sum(T) - complete


def configuration_rows(graph, H, C, T):
    rows = []
    for ray, label in zip(graph.rays, graph.label):
        if label[0] == 'hex':
            i, h = label[1:]
            take = H[i] == 2 or (H[i] == 1 and h in [(1, 0), (-1, 0)])
        else:
            c, t = label[1:]
            take = T[c] if t == 0 else C[c]
        if take:
            rows.append(ray)
    return rows


def actual_rank(graph, H, C, T):
    rows = configuration_rows(graph, H, C, T)
    return G.rref(rows)[0] if rows else 0


def check_lifted_automorphism(graph, pi, reflected=False):
    """Map each actual ray, with signs on reversed hexagon factors."""
    # This also handles j=2: its two parallel edges are swapped by rotation.
    j = graph.j
    if reflected:
        hp = [(-i-1) % j for i in range(j)]
        signs = [-1]*j
    else:
        hp = [(i+1) % j for i in range(j)]
        signs = [1]*j
    rayset = set(graph.rays)
    for ray, label in zip(graph.rays, graph.label):
        out = [0]*graph.n
        for i in range(j):
            out[2*hp[i]] = signs[i]*ray[2*i]
            out[2*hp[i]+1] = signs[i]*ray[2*i+1]
        for c in range(j):
            for t in range(2):
                out[graph.off[pi[c]]+t] = ray[graph.off[c]+t]
        assert tuple(out) in rayset
        if label[0] == 'base' and label[2] == 0:
            target = graph.rays[graph.label.index(('base', pi[label[1]], 0))]
            assert tuple(out) == target
    # Each block is sent to a block, hexagon adjacency survives +/- identity,
    # and each base omitted-ray choice is permuted. Thus maximal cones survive.


def rank_gates():
    checked = 0
    for graph in [G.path(1), G.path(2), G.star(3, 3), cycle(2), cycle(3)]:
        # Complete configurations on path1/2; deterministic sample on larger graphs.
        for i, flat in enumerate(product(range(3), repeat=graph.j)):
            for z, bits in enumerate(product(range(2), repeat=2*graph.k)):
                if graph.j >= 3 and (i*137 + z*31) % 97:
                    continue
                C, T = bits[:graph.k], bits[graph.k:]
                assert graph_rank(graph, flat, C, T) == actual_rank(graph, flat, C, T)
                checked += 1
    print(f'General graph rank formula: {checked} rational-rank comparisons passed.', flush=True)
    return checked


def cycle_gates():
    records = []
    for j in range(1, 41):
        Z, a, b, c = cycle_weights(j)
        s = [Z/2 - Z*a, Z*b - Z/4, Z*c-F(2, 3)*Z, F(5, 6)*Z-Z*c]
        if j >= 2:
            assert a < F(1, 2) and b > F(1, 4) and F(2, 3) < c < F(5, 6)
        records.append({'j': j, 'Z': str(Z), 'a': str(a), 'b': str(b), 'c': str(c),
                        'positive_sequences': list(map(str, s))})
        if j <= 3:
            print(f'j={j} exact sign sequences: {s}', flush=True)
    seeds = [[F(x) for x in row['positive_sequences']] for row in records[:3]]
    assert all(x >= 0 for row in seeds for x in row)
    assert all(x > 0 for row in seeds[1:] for x in row)
    assert all(x > 12*y for x, y in zip(seeds[2], seeds[1]))
    assert F(27, 2)-F(2491, 240)/12 > 12
    for i in range(3, len(records)):
        rows = [[F(x) for x in row['positive_sequences']] for row in records[i-3:i+1]]
        for q in range(4):
            assert rows[3][q] == F(27, 2)*rows[2][q] - F(2491, 240)*rows[1][q] + F(127, 864)*rows[0][q]
    for j in range(2, 9):
        graph = cycle(j)
        check_lifted_automorphism(graph, [(c+1) % j for c in range(j)])
        check_lifted_automorphism(graph, [(-c) % j for c in range(j)], reflected=True)
        Z, a, b, c = cycle_weights(j)
        for h, C, T in product(range(3), range(2), range(2)):
            Hs, Cs, Ts = [h]*j, [C]*j, [T]*j
            rank = graph_rank(graph, Hs, Cs, Ts)
            assert rank == actual_rank(graph, Hs, Cs, Ts)
            weight = j*((0, 2*a, 2*a+4*b)[h] + (2*C+T)*c)
            if 0 < rank < graph.n:
                assert rank > weight
        if j <= 3:
            deg, D = graph.degrees()
            for d, label in zip(deg, graph.label):
                expected = c if label[0] == 'base' else (a if label[2] in [(1, 0), (-1, 0)] else b)
                assert graph.n*d/D == expected
            assert D == __import__('math').factorial(graph.n)*Z
    print('Cycle gates: trace degrees match independent polynomial integration for j=2,3;', flush=True)
    print('actual lifted-ray symmetries and orbit-union ranks checked for j=2..8;', flush=True)
    print('four trace recurrences and inequalities checked exactly through j=40.', flush=True)
    return records


def tree_weights(graph):
    """Exact polynomial messages; independently checked against old certificates."""
    adj = [[] for _ in graph.b]
    for e, (p, m) in enumerate(graph.edges):
        adj[p].append((m, e, 1)); adj[m].append((p, e, -1))

    def moment(poly, k):
        return sum(a*P.I(i+k) for i, a in enumerate(poly))

    def factor(c, parent, power):
        children = [(v, e, sign) for v, e, sign in adj[c] if v != parent]
        boundary_sign = next((sign for v, e, sign in adj[c] if v == parent), 0)
        polys = [message(v, c) for v, _, _ in children]
        out = []
        for q in range(power+1 if parent is not None else 1):
            val = F(0)
            for ks in product(range(power-q+1), repeat=len(children)):
                rem = power-q-sum(ks)
                if rem < 0: continue
                term = F((graph.b[c]+1)**rem, math.factorial(rem)*math.factorial(q))
                term *= boundary_sign**q
                for k, poly, (_, _, sign) in zip(ks, polys, children):
                    term *= F(sign**k, math.factorial(k))*moment(poly, k)
                val += term
            out.append(val)
        return tuple(out)

    @lru_cache(None)
    def message(c, parent):
        return factor(c, parent, graph.b[c])

    Z = factor(0, None, graph.b[0])[0]
    base = [factor(c, None, graph.b[c]-1)[0]/Z for c in range(graph.k)]
    hexa = []
    for p, m in graph.edges:
        u, v = message(p, m), message(m, p)
        vals = {kind: sum(x*y*fn(i+k) for i, x in enumerate(u) for k, y in enumerate(v))/Z
                for kind, fn in P.EDGE.items()}
        hexa.append(vals)
    delta = [base[l[1]] if l[0] == 'base' else hexa[l[1]][G.KIND[l[2]]] for l in graph.label]
    assert sum(delta) == graph.n
    return Z, base, hexa, delta


def tree_deficit(graph, base, hexa):
    """Minimum over all locally invariant configurations, by exact tree recursion.

    States (a,e,f): open component entirely C=T=1; nonempty selection;
    all processed nontrivial representations needed for full rank selected.
    """
    adj = [[] for _ in graph.b]
    for e, (p, m) in enumerate(graph.edges):
        adj[p].append((m, e)); adj[m].append((p, e))

    def update(dic, key, value):
        if key not in dic or value < dic[key]: dic[key] = value

    def visit(c, parent):
        table = {}
        for C, T in product(range(2), repeat=2):
            state = (bool(C*T), bool(C or T), bool(C or (graph.b[c] == 1 and T)))
            update(table, state, (graph.b[c]*C+T)*(1-base[c]))
        for child, edge in adj[c]:
            if child == parent: continue
            sub = visit(child, c); nxt = {}
            h = hexa[edge]
            weights = [0, h['r']+h['-r'], h['r']+h['-r']+2*h['neg']+2*h['pos']]
            for (a,e,f), cost in table.items():
                for (aa,ee,ff), ccost in sub.items():
                    for H in range(3):
                        state = (a and aa if H == 0 else a, e or ee or H>0, f and ff and H==2)
                        update(nxt, state, cost+ccost+H-weights[H]-(int(aa) if H>0 else 0))
            table = nxt
        return table
    table = visit(0, None)
    return min(value-int(a) for (a,e,f), value in table.items() if e and not f)


def spider_gates():
    result = []
    cases = [(1,1,1),(1,1,2),(1,2,2),(1,1,3),(2,2,2),(2,2,3),(2,3,5)]
    cases += [(L,L,L) for L in [3,4,8,16,32]]
    for legs in cases:
        graph = G.spider(*legs)
        Z, base, hexa, delta = tree_weights(graph)
        deficit = tree_deficit(graph, base, hexa)
        if legs in [(1,1,1),(1,1,2),(1,2,2),(1,1,3)]:
            degrees, D = graph.degrees()
            assert delta == [graph.n*d/D for d in degrees]
        if legs == (2,2,2):
            cert = json.loads((Path(__file__).resolve().parents[1]/'certificates/spider222.json').read_text())
            assert graph.rays == list(map(tuple, cert['rays']))
            assert delta == list(map(F, cert['delta']))
        row = {'legs':legs, 'n':graph.n, 'minimum_invariant_deficit':str(deficit),
               'deficit_decimal':float(deficit), 'base_weights':list(map(str,base))}
        result.append(row)
        print(f'Spider {legs}: n={graph.n}, minimum invariant deficit={float(deficit):+.12g}', flush=True)
    assert result[0]['deficit_decimal'] > 0
    assert all(r['deficit_decimal'] < 0 for r in result[1:4])
    assert all(r['deficit_decimal'] > 0 for r in result[4:])
    return result


def tree_dp_gate():
    checked = 0
    for graph in [G.path(1), G.path(2), G.star(3, 3)]:
        Z, base, hexa, delta = tree_weights(graph)
        answer = tree_deficit(graph, base, hexa)
        brute = None
        for H in product(range(3), repeat=graph.j):
            weight_h = sum((0, h['r']+h['-r'], h['r']+h['-r']+2*h['neg']+2*h['pos'])[x]
                           for x, h in zip(H, hexa))
            for bits in product(range(2), repeat=2*graph.k):
                C, T = bits[:graph.k], bits[graph.k:]
                rank = graph_rank(graph, H, C, T)
                if not 0 < rank < graph.n: continue
                margin = rank-weight_h-sum((b*c+t)*w for b,c,t,w in zip(graph.b,C,T,base))
                if brute is None or margin < brute: brute = margin
                checked += 1
        assert brute == answer
    print(f'Tree minimum-deficit recurrence: exhaustive comparison on {checked} proper configurations passed.', flush=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    count = rank_gates()
    records = cycle_gates()
    (OUT/'cycles.json').write_text(json.dumps({'rank_comparisons': count, 'cycles': records}, indent=2)+'\n')
    spiders = spider_gates()
    (OUT/'spiders.json').write_text(json.dumps(spiders, indent=2)+'\n')
    tree_dp_gate()
    print('All bounded exact gates passed. All-length implication is the written induction.', flush=True)


if __name__ == '__main__':
    main()
