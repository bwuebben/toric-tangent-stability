"""Exact constants for the two-junction long-bridge argument.

The infinite implication is proved in two_junctions_proof.md. This program
checks the bulk inequalities and rational error budget, not a length census.
Uses only the standard library and unchanged reviewed transfer modules.
"""
from fractions import Fraction as F
from itertools import product
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import path_ends as PE
import balanced_spider_proof as BS
import unequal_spider_proof as US

OUT = Path(__file__).resolve().parents[1] / 'output/research_two_junction_uniform'


def constants(q=16, P=40):
    B, D, tau = F(81), F(11757, 10000), F(2, 7)
    x, z = D*tau**(q-1), D*tau**(2*q-1)
    assert 2*x < 1 and z < 1
    ends = B*(16+36/(1-tau))*(z+z*z)
    bridge = B*(18*q*(x+x*x)+18*x/(1-tau)*(1+2*x))
    error = ends+bridge
    y = 2*D*tau**(P-1)
    ray_error = B*(y+y*y)
    R = BS.right_chains(P)[P]
    L = US.left_chain(PE.end_vec(1, -1), P)[P]
    w = PE.hex_weights(L, R)[0]
    assert w[2]-6*ray_error > F(187, 100)
    assert w[2]+6*ray_error < F(15, 8)
    assert w[1]+2*ray_error < F(14, 25)
    # The true bulk identity is exact; impose it for the symbolic check.
    d, t = (4-w[2])/3, w[2]/2
    transfers = []
    for a, C, T, H in product((0, 1), (0, 1), (0, 1), range(3)):
        anew = C*T*(a if H == 0 else 1)
        value = (2*C+T)*(1-d)+H-w[H]-a*(H > 0)+t*(a-anew)
        zero = (a, C, T, H) in ((0, 0, 0, 0), (1, 1, 1, 2))
        assert value == 0 if zero else value > F(1, 16)
        transfers.append(dict(a=a, C=C, T=T, H=H, anew=anew, cost=value, zero=zero))
    joins = []
    for a, b, H in product((0, 1), (0, 1), range(3)):
        value = t*(a+b)-a*b if H == 0 else H-w[H]+(t-1)*(a+b)
        zero = (a, b, H) in ((0, 0, 0), (1, 1, 2))
        assert value == 0 if zero else value > F(1, 16)
        joins.append(dict(a=a, b=b, H=H, cost=value, zero=zero))
    mu = F(2698, 100000)
    lower = min(mu, F(1, 16))-2*error
    assert error < F(206, 10**6)
    assert lower > F(2656, 100000)
    return dict(q=q, P=P, minimum_bridge=2*q+1, outer_minimum=2,
                B=B, D=D, tau=tau, x=x, z=z, end_error=ends,
                bridge_error=bridge, total_error=error,
                bulk_proxy=w, bulk_ray_error=ray_error,
                inherited_margin=mu, bulk_gap=F(1, 16),
                full_endpoint_error=error/2, uniform_lower=lower,
                transfers=transfers, joins=joins)


if __name__ == '__main__':
    data = constants()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT/'constants.json').write_text(json.dumps(BS.serialize(data), indent=2)+'\n')
    print('Bulk: 24 transfers and 12 joins; exactly two zero choices each.')
    print('All remaining bulk costs exceed 1/16 at the rational proxy.')
    print('Certified true-weight bounds: 1.87 < w2 < 1.875; w1 < 0.56.')
    print('End error:', float(data['end_error']))
    print('Bridge error:', float(data['bridge_error']))
    print('Uniform error:', float(data['total_error']))
    print('Final lower bound:', float(data['uniform_lower']))
    print('Exact constants PASS; theorem requires the written gluing argument and review.')
