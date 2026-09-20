"""Reproduce explicit sufficient subdivision thresholds in the G3 proof.

No stability census is performed. The theorem is the induction in
research_notes/Toric interactions execution/subcubic_gluing_proof.md.
"""
from fractions import Fraction as F
from pathlib import Path
import json

MU0, B, D, TAU = F(2698, 100000), F(81), F(11757, 10000), F(2, 7)


def geometric_error(junctions, q):
    assert junctions >= 2 and q >= 1
    x = D*TAU**(q-1)
    return B*F(156*junctions+272, 5)*x*(1+2*x)


def threshold(junctions):
    assert junctions >= 1
    margin = MU0*(F(1, 2)+F(1, 2*junctions))
    if junctions == 1:
        return dict(junctions=1, length=2, margin=margin)
    q = 1
    budget = MU0/(4*junctions*(junctions-1))
    while 2*D*TAU**(q-1) > 1 or geometric_error(junctions, q) > budget:
        q += 1
    error = geometric_error(junctions, q)
    previous_margin = MU0*(F(1, 2)+F(1, 2*(junctions-1)))
    assert previous_margin-2*budget == margin
    assert previous_margin-2*error >= margin
    return dict(junctions=junctions, q=q, length=2*q+1,
                margin=margin, geometric_error=error, error_budget=budget)


def serialize(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: serialize(v) for k, v in value.items()}
    if isinstance(value, list):
        return [serialize(v) for v in value]
    return value


if __name__ == '__main__':
    rows = [threshold(n) for n in (1, 2, 3, 4, 5, 10, 25, 50, 100, 1000)]
    out = Path(__file__).resolve().parents[1]/'output/subcubic_thresholds'
    out.mkdir(parents=True, exist_ok=True)
    (out/'thresholds.json').write_text(json.dumps(serialize(rows), indent=2)+'\n')
    for row in rows:
        print('junctions', row['junctions'], 'minimum edge length', row['length'],
              'invariant margin >=', float(row['margin']))
    print('Exact threshold arithmetic PASS; the all-tree assertion is the written induction.')
