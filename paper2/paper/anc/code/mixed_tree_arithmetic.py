#!/usr/bin/env python3
"""Exact finite arithmetic printed in the mixed-blowup tree proofs.

Uses the split-at-one coefficient identity in the manuscript. No enumeration
of tree shapes or extrapolation in path length is part of this computation.
Run with --check-json FILE to compare a saved record outside the checkout.
"""
import argparse
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import json


def multiply(p, q, degree):
    result = [F(0)]*(degree+1)
    for i,a in enumerate(p):
        for j,b in enumerate(q):
            if i+j <= degree:
                result[i+j] += a*b
    return tuple(result)


def Q(p, degree):
    return (sum(p),) + tuple(sum((a*F(factorial(k),factorial(k+j+1))
          for k,a in enumerate(p)),F(0)) for j in range(degree))


def series(inputs, degree, offset=1):
    s = tuple(F(offset)**j/factorial(j) for j in range(degree+1))
    for p in inputs:
        s = multiply(s,Q(p,degree),degree)
    return s


def outgoing(inputs):
    d = len(inputs)+1
    s = series(inputs,d)
    return tuple(s[d-j]/factorial(j) for j in range(d+1))


def root_record(inputs, degree):
    s = series(inputs,degree)
    return {'volume':s[degree], 'facet':s[degree-1],
            'degree':s[degree-1]/s[degree], 'deficit':1-s[degree-1]/s[degree]}


def compute():
    ell = (F(1),F(1))
    p = tuple(map(F,(8,12,6,1)))
    closure = []
    for d in (2,3):
        g = outgoing([p]*(d-1))
        a = g+(F(0),)*(4-len(g))
        numerator = (2*a[1]-3*a[0],4*a[2]-2*a[1],6*a[3]-a[2])
        assert all(x<0 for x in numerator)
        closure.append({'degree':d, 'polynomial':g, 'derivative_numerator':numerator})
    subcubic = root_record([p]*3,3)
    assert subcubic['degree']==F(25950132,26414113)
    gs = {d:outgoing([ell]*(d-1)) for d in (3,4,5)}
    h = tuple(120*x for x in gs[4])
    dh = tuple(j*h[j] for j in range(1,len(h)))
    ddh = tuple(j*dh[j] for j in range(1,len(dh)))
    lc = tuple(a-b for a,b in zip(multiply(dh,dh,6),multiply(h,ddh,6)))
    assert all(x>0 for x in lc)
    degree_four = root_record([gs[4]]*4,4)
    assert degree_four['degree']==F(40090063176909393024,39498346443818421703)
    phi2 = (F(2),F(2),F(1,2))
    degree_five = root_record([phi2]*4+[ell],5)
    assert degree_five['degree']==F(6244655501,6233518181)
    spider = root_record([outgoing([ell])]*5,5)
    assert spider['deficit']==-F(304924770977,12705819178443)
    signs, endpoints, pairs = [], [], []
    for d in (3,4,5):
        values=[]
        for k in range(6):
            monomial = (F(0),)*k+(F(1),)
            s = series([ell]*(d-1)+[outgoing([monomial])],d)
            values.append(s[d]-s[d-1])
        assert all(x>0 for x in values)
        signs.append({'degree':d,'values':values})
        gap=sum((1-j)*a for j,a in enumerate(gs[d]))
        assert gap>0
        endpoints.append({'degree':d,'polynomial':gs[d], 'value_minus_derivative':gap})
        for e in range(d,6):
            left=root_record([ell]*(d-1)+[gs[e]],d)
            right=root_record([ell]*(e-1)+[gs[d]],e)
            assert left['volume']==right['volume']
            gap=min(left['deficit'],right['deficit'])
            assert gap>0
            pairs.append({'degrees':[d,e],'minimum_deficit':gap})
    return {'scope':'printed finite arithmetic; uniform conclusions use the manuscript proofs',
            'reference':{'polynomial':p,'Q_coefficients':Q(p,3),'closure':closure,'root':subcubic},
            'degree_four':{'polynomial':gs[4],'log_concavity_numerator':lc,
                           'Q_coefficients':Q(gs[4],4),'root':degree_four},
            'degree_five':degree_five, 'five_arms_length_two':spider,
            'two_centres':{'endpoints':endpoints,'functional_signs':signs,'adjacent_pairs':pairs}}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-json',type=Path)
    parser.add_argument('--write-json',type=Path)
    args=parser.parse_args()
    record=json.loads(json.dumps(compute(),default=str))
    if args.check_json:
        assert record==json.loads(args.check_json.read_text()), 'Recorded arithmetic differs'
    if args.write_json:
        args.write_json.parent.mkdir(parents=True,exist_ok=True)
        args.write_json.write_text(json.dumps(record,indent=2)+'\n')
    print('PASS: subcubic closure and root; degree-four and degree-five obstructions; '
          'two-centre endpoint signs, 18 functional signs, and six adjacent pairs.')


if __name__=='__main__':
    main()
