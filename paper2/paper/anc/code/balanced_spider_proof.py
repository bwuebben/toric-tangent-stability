"""Finite rational data for an all-length balanced trivalent-spider proof.

This script implements the printed one-arm/junction recurrence and comparisons.
It does not modify the existing graph or path proof scripts.
"""
from fractions import Fraction as F
from itertools import product, permutations
from pathlib import Path
import json
import math
import random
import sys
import time

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import path_transfer as PT
import path_ends as PE
import path_all_j as PA
import research_graph_targets as RG
import tree_family as TF

OUT = Path(__file__).resolve().parents[1]/'output/research_balanced_spiders'
BOOL = list(product((False, True), repeat=2))


def normalize(v):
    assert v[0] > 0
    return [x/v[0] for x in v]


def right_chains(top):
    rows = [PE.end_vec(1, 1)]
    for _ in range(top): rows.append(normalize(PE.mv(PE.TR, rows[-1])))
    return rows


def junction_polynomial(h, power=3):
    """Integrate (4-a-b-c)^power/power! against h(b)h(c)dnu(b)dnu(c)."""
    moments = [sum(x*PT.I(i+k) for i,x in enumerate(h)) for k in range(4)]
    out = []
    for a in range(4):
        val = F(0)
        for b in range(power-a+1):
            for c in range(power-a-b+1):
                d = power-a-b-c
                val += F((-1)**(a+b+c)*4**d,
                         math.factorial(a)*math.factorial(b)*math.factorial(c)*math.factorial(d))*moments[b]*moments[c]
        out.append(val)
    return out


def data_from_chains(length, R):
    h = R[length-1]
    J, Jdrop = junction_polynomial(h), junction_polynomial(h,2)
    L = [normalize(J)]
    for _ in range(length-1): L.append(normalize(PE.vm(L[-1],PE.T)))
    central = PE.dot(PE.vm(Jdrop,PE.W),h)/PE.dot(PE.vm(J,PE.W),h)
    hw = [PE.hex_weights(L[i],R[length-1-i])[0] for i in range(length)]
    bw = [None]+[PE.base_weight(L[c-1],R[length-1-c]) for c in range(1,length)]
    leaf = PE.dot(PE.vm(L[-1],PE.W),PE.end_vec(1,1,1))/PE.dot(PE.vm(L[-1],PE.W),PE.end_vec(1,1))
    assert 3*(sum(w[2] for w in hw)+3*sum(bw[1:])+2*leaf)+4*central == 12*length
    return hw,bw,leaf,central


def put(table,key,value):
    if key not in table or value < table[key]: table[key]=value


def arm_states(hw,bw,leaf):
    """Read leaf to junction; output G/g/u, nonempty and full flags."""
    length=len(hw); table={}
    for C,T in BOOL:
        base=F(int(C)+int(T and not C))-(int(C)+int(T))*leaf
        for H in range(3):
            link=C and T
            cost=base+H-hw[-1][H]
            if H>0: state='G'
            elif link: cost+=1; state='g'
            else: state='u'
            put(table,(state,bool(C or T or H),bool((C or T) and H==2)),cost)
    for position in range(length-2,-1,-1):
        nxt={}; base=bw[position+1]
        for (state,e,f),cost in table.items():
            for C,T in BOOL:
                link=C and T
                for H in range(3):
                    v=cost+2*int(C)+int(T and not C)-(2*int(C)+int(T))*base+H-hw[position][H]
                    if H>0: v+=int(link and state=='u'); ns='G'
                    elif link: v+=1; ns='g' if state in 'Gg' else 'u'
                    else: ns='u'
                    put(nxt,(ns,bool(e or C or T or H),bool(f and C and H==2)),v)
        table=nxt
    return table


def junction_table(hw,bw,leaf,central):
    arms=arm_states(hw,bw,leaf); rows=[]
    for (state,e,f),cost in sorted(arms.items()):
        for C,T in BOOL:
            empty=not(e or C or T); full=f and C
            rankadd=3*int(C)+int(T and not C)+int(C and T and state=='u')
            deficit=3*cost+rankadd-(3*int(C)+int(T))*central
            rows.append(dict(state=state,nonempty=e,arm_full=f,C=C,T=T,
                             rankadd=rankadd,excluded='empty' if empty else 'full' if full else None,
                             deficit=deficit,arm_cost=cost))
    return rows


def minimum(data):
    rows=junction_table(*data)
    return min(r['deficit'] for r in rows if not r['excluded'])


def model_data(m,q,P,R):
    length=2*q+m; J=junction_polynomial(R[P]); Jdrop=junction_polynomial(R[P],2)
    L=[normalize(J)]
    for _ in range(max(P,length)): L.append(normalize(PE.vm(L[-1],PE.T)))
    def indices(position):
        if position<q:return 'junction'
        if position<q+m:return 'bulk'
        return 'leaf'
    hw=[];bw=[None]
    for i in range(length):
        region=indices(i)
        hw.append(PE.hex_weights(L[i] if region=='junction' else L[P],
                                 R[length-1-i] if region=='leaf' else R[P])[0])
    for c in range(1,length):
        # The rank recurrence runs from leaf to junction: one period is the
        # base at c followed by hexagon c-1. Zone membership must follow that
        # common period index so the central block has exactly m bulk periods.
        region=indices(c-1)
        bw.append(PE.base_weight(L[c-1] if region=='junction' else L[P],
                                 R[length-1-c] if region=='leaf' else R[P]))
    leaf=PE.dot(PE.vm(L[P],PE.W),PE.end_vec(1,1,1))/PE.dot(PE.vm(L[P],PE.W),PE.end_vec(1,1))
    central=PE.dot(PE.vm(Jdrop,PE.W),R[P])/PE.dot(PE.vm(J,PE.W),R[P])
    return hw,bw,leaf,central


def errors(q=20,P=40):
    B=F(81);D=F(11757,10000);tau=F(2,7)
    assert sum(D**r/math.factorial(r) for r in range(8))>F(81,25)
    x=D*tau**(q-2);y=3*D*tau**(P-1);z=3*D*tau**(2*q-2)
    s=2*x/(1-tau)
    zoning=3*B*(9*s*(1+2*x)+2*(9*q+4)*(x+x*x))
    junction=B*(54*q+1)*(z+z*z)
    proxy=B*(27*(2*q+3)+1)*(y+y*y)
    return dict(zoning=zoning,junction=junction,proxy=proxy,total=zoning+junction)


def gates():
    # Direct fan action under arbitrary arm permutations, on every actual ray.
    for length in [1,2,3]:
        graph=TF.spider(length,length,length); lookup=set(graph.rays)
        for pi in permutations(range(3)):
            for ray in graph.rays:
                out=list(ray)
                for a in range(3):
                    for i in range(length):
                        edge=a*length+i; image=pi[a]*length+i
                        out[2*image:2*image+2]=ray[2*edge:2*edge+2]
                        c=1+a*length+i; v=1+pi[a]*length+i
                        for t in range(graph.b[c]):out[graph.off[v]+t]=ray[graph.off[c]+t]
                assert tuple(out) in lookup
        # Exact degrees independently by the previously tested tree messages.
        _,bases,hexes,_=RG.tree_weights(graph)
        data=data_from_chains(length,right_chains(length));hw,bw,leaf,central=data
        assert central==bases[0] and leaf==bases[length]
        for i in range(length):
            h=hexes[i]
            assert hw[i]==(0,h['r']+h['-r'],h['r']+h['-r']+2*h['neg']+2*h['pos'])
            if i:assert bw[i]==bases[i]
        # Exhaustive symmetric configurations for L<=2; deterministic L3 sample.
        rng=random.Random(932+length); expected=None; count=0
        configurations=product(range(3),repeat=length)
        for hs in configurations:
            for leafC,leafT in BOOL:
                for mids in product(BOOL,repeat=length-1):
                    for centerC,centerT in BOOL:
                        if length==3 and rng.randrange(13):continue
                        Cs=[int(centerC)];Ts=[int(centerT)];Hs=[]
                        for a in range(3):
                            Hs.extend(hs)
                            Cs.extend([int(x[0]) for x in mids]+[int(leafC)])
                            Ts.extend([int(x[1]) for x in mids]+[int(leafT)])
                        actual=RG.actual_rank(graph,Hs,Cs,Ts)
                        # Independently scan this one fixed arm with integer ranks.
                        rank=int(leafC)+int(leafT and not leafC);state=None
                        for i in range(length-1,-1,-1):
                            if i==length-1:C,T=leafC,leafT
                            else:
                                C,T=mids[i];rank+=2*int(C)+int(T and not C)
                            link=C and T;H=hs[i]
                            if H:rank+=H+int(link and state=='u');state='G'
                            elif link:rank+=1;state='g' if state in [None,'G','g'] else 'u'
                            else:state='u'
                        proposed=3*rank+3*int(centerC)+int(centerT and not centerC)+int(centerC and centerT and state=='u')
                        assert actual==proposed
                        if 0<actual<graph.n:
                            weight=3*(sum(hw[i][hs[i]] for i in range(length))+sum(bw[i+1]*(2*C+T) for i,(C,T) in enumerate(mids))+leaf*(leafC+leafT))+central*(3*centerC+centerT)
                            margin=actual-weight
                            if expected is None or margin<expected:expected=margin
                        count+=1
        if length<=2:assert expected==minimum(data)
        print(f'Gate L={length}: fan permutations, degrees, {count} direct symmetric ranks passed.',flush=True)


def serialize(x):
    if isinstance(x,F):return str(x)
    if isinstance(x,dict):return {k:serialize(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):return [serialize(v) for v in x]
    return x


def main():
    started=time.monotonic();OUT.mkdir(parents=True,exist_ok=True);gates()
    q,P=20,40;R=right_chains(max(P,2*q+3));record=dict(q=q,proxy=P,finite=[],models=[])
    for length in range(2,2*q+4):
        value=minimum(data_from_chains(length,R));assert value>0
        record['finite'].append(dict(length=length,minimum=value))
    worst=min(record['finite'],key=lambda row:row['minimum'])
    print(f'Finite L=2..{2*q+3}: minimum {float(worst["minimum"]):.12g} at L={worst["length"]}',flush=True)
    full=[]
    for m in range(4):
        data=model_data(m,q,P,R);value=minimum(data);hw,bw,leaf,central=data;length=len(hw)
        fullcost=12*length-(3*(sum(h[2] for h in hw)+3*sum(bw[1:])+2*leaf)+4*central)
        full.append(abs(fullcost));rows=junction_table(*data)
        record['models'].append(dict(m=m,minimum=value,fullcost=fullcost,junction_table=rows))
        print(f'Model m={m}: min {float(value):.12g}; full cost {float(fullcost):+.3e}; junction rays {float(central):.12g}',flush=True)
        if m==0:
            for row in rows:
                print(f'  junction {row["state"]}, e={int(row["nonempty"])}, f={int(row["arm_full"])}, C={int(row["C"])}, T={int(row["T"])}: {float(row["deficit"]):+.12f} {row["excluded"] or "proper"}',flush=True)
    Lproxy=normalize(PE.end_vec(1,-1))
    for _ in range(P):Lproxy=normalize(PE.vm(Lproxy,PE.T))
    hw=PE.hex_weights(Lproxy,R[P])[0];bw=PE.base_weight(Lproxy,R[P])
    (gamma,_,_),trans=PA.gamma_cyc(PA.period_automaton(hw,bw))
    errs=errors(q,P);record['errors']=errs;record['cycle_gap']=gamma
    record['transitions']={str(k):v[0] for k,v in trans.items()}
    full_bound=max(full)+errs['proxy']
    lower=min(min(x['minimum'] for x in record['models'])-errs['proxy'],3*gamma-errs['proxy']-full_bound)-errs['total']
    record['full_bound']=full_bound;record['all_length_lower_bound']=lower
    record['junction_lower_bounds']=[
        dict(state=row['state'],nonempty=row['nonempty'],arm_full=row['arm_full'],
             C=row['C'],T=row['T'],excluded=row['excluded'],
             lower_six_decimals=F(0))
        for i,row in enumerate(record['models'][0]['junction_table'])]
    for i,row in enumerate(record['junction_lower_bounds']):
        lo=min(mod['junction_table'][i]['deficit'] for mod in record['models'])
        row['lower_six_decimals']=F((lo*10**6).__floor__(),10**6)
    print('Error bounds:',{k:float(v) for k,v in errs.items()},flush=True)
    print(f'Positive arm-cycle gap {float(gamma):.12g}; global lower bound L>=40: {float(lower):.12g}',flush=True)
    assert lower>F(8849,100000)
    assert worst['minimum']>F(876,10000)
    assert all(mod['minimum']>F(88516,10**6) for mod in record['models'])
    assert errs['total']<F(181,10**7) and errs['proxy']<F(21,10**17)
    assert gamma>F(1259,10000) and max(full)<F(2,10**24)
    printed=["0.290258","0.514907","0.964204","0.188853",
             "0.101404","0.326053",None,None,
             "2.767585","2.992234","3.441531","2.666180",
             None,"0.224648","0.673946","0.898595",
             "0.088516","0.313165","0.762462","0.987111"]
    for mod in record['models']:
        for row,bound in zip(mod['junction_table'],printed):
            if bound is not None:assert row['deficit']>=F(bound)
    (OUT/'certificate.json').write_text(json.dumps(serialize(record),indent=2)+'\n')
    print(f'All finite inputs and comparisons passed in {time.monotonic()-started:.2f}s.',flush=True)


if __name__=='__main__':main()
