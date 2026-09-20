"""Exact asymmetric arm/junction calculations; research, not yet a theorem.

The initial diagnostic is declared in research_notes/Toric interactions
execution/unequal_arms.md. Reviewed balanced/path modules are imported unchanged.
"""
from fractions import Fraction as F
from itertools import product, combinations_with_replacement
from pathlib import Path
import json
import math
import sys
import time
from functools import lru_cache

if hasattr(sys, 'set_int_max_str_digits'):
    sys.set_int_max_str_digits(0)
sys.path.insert(0, str(Path(__file__).resolve().parent))
import balanced_spider_proof as BS
import research_graph_targets as RG
import tree_family as TF
import path_ends as PE
import path_transfer as PT
import path_all_j as PA

OUT = Path(__file__).resolve().parents[1] / 'output/research_unequal_spiders'
BOOL = tuple(product((False, True), repeat=2))


def mixed_junction(h, k, power=3):
    mh = [sum(x*PT.I(i+t) for i,x in enumerate(h)) for t in range(4)]
    mk = [sum(x*PT.I(i+t) for i,x in enumerate(k)) for t in range(4)]
    out = [F(0)]*4
    for a in range(power+1):
        for b in range(power-a+1):
            for c in range(power-a-b+1):
                d = power-a-b-c
                out[a] += F((-1)**(a+b+c)*4**d,
                            math.factorial(a)*math.factorial(b)*math.factorial(c)*math.factorial(d))*mh[b]*mk[c]
    return out


def left_chain(J, top):
    L = [BS.normalize(J)]
    for _ in range(top):
        L.append(BS.normalize(PE.vm(L[-1], PE.T)))
    return L


def exact_data(lengths, R):
    hs = [R[n-1] for n in lengths]
    arms = []
    centers = []
    for a,n in enumerate(lengths):
        others = [hs[b] for b in range(3) if b != a]
        J = mixed_junction(*others)
        Jminus = mixed_junction(*others, power=2)
        centers.append(PE.dot(PE.vm(Jminus,PE.W),hs[a])/PE.dot(PE.vm(J,PE.W),hs[a]))
        L = left_chain(J,n-1)
        hw = [PE.hex_weights(L[i],R[n-1-i])[0] for i in range(n)]
        bw = [None]+[PE.base_weight(L[c-1],R[n-1-c]) for c in range(1,n)]
        leaf = PE.dot(PE.vm(L[-1],PE.W),PE.end_vec(1,1,1))/PE.dot(PE.vm(L[-1],PE.W),PE.end_vec(1,1))
        arms.append((hw,bw,leaf))
    assert centers[0] == centers[1] == centers[2]
    central = centers[0]
    assert sum(sum(h[2] for h in hw)+3*sum(bw[1:])+2*leaf for hw,bw,leaf in arms)+4*central == 4*sum(lengths)
    return arms,central


def central_rows(arms, central):
    tables = [BS.arm_states(*arm) for arm in arms]
    rows = []
    for states in product(*(sorted(table) for table in tables)):
        armcost = sum(table[state] for table,state in zip(tables,states))
        for C,T in BOOL:
            empty = not(C or T or any(s[1] for s in states))
            full = C and all(s[2] for s in states)
            rankadd = 3*int(C)+int(T and not C)+int(C and T and any(s[0]=='u' for s in states))
            rows.append(dict(states=states,C=C,T=T,excluded='empty' if empty else 'full' if full else None,
                             deficit=armcost+rankadd-(3*int(C)+int(T))*central))
    return rows


def minimum(data):
    return min((row for row in central_rows(*data) if not row['excluded']), key=lambda row:row['deficit'])


def check_against_tree(lengths,data):
    graph = TF.spider(*lengths)
    _,bases,hexes,_ = RG.tree_weights(graph)
    arms,central = data
    assert bases[0] == central
    pos = 0
    for n,(hw,bw,leaf) in zip(lengths,arms):
        assert bases[pos+n] == leaf
        for i in range(n):
            h=hexes[pos+i]
            assert hw[i] == (0,h['r']+h['-r'],h['r']+h['-r']+2*h['neg']+2*h['pos'])
            if i:
                assert bw[i] == bases[pos+i]
        pos += n
    assert minimum(data)['deficit'] == RG.tree_deficit(graph,bases,hexes)


def diagnostic():
    OUT.mkdir(parents=True,exist_ok=True)
    cases = list(combinations_with_replacement(range(2,7),3))
    cases += [(2,2,12),(2,3,12),(2,6,12),(3,3,12),(6,6,12),
              (2,12,12),(3,12,12),(6,12,12),(12,12,12),(2,2,40)]
    controls = [(1,1,1),(1,1,2),(1,2,2),(1,1,3)]
    R = BS.right_chains(40)
    rows=[]
    for lengths in controls+cases:
        data=exact_data(lengths,R)
        check_against_tree(lengths,data)
        row=minimum(data)
        row.update(lengths=lengths,central=data[1],kind='control' if lengths in controls else 'target')
        rows.append(row)
        print(lengths, float(row['deficit']),row['states'],row['C'],row['T'],flush=True)
    record=dict(status='diagnostic only',rows=rows,finite_target_count=len(cases))
    (OUT/'diagnostic.json').write_text(json.dumps(BS.serialize(record),indent=2)+'\n')
    print('Worst target:',min((r for r in rows if r['kind']=='target'),key=lambda r:r['deficit']),flush=True)


def error_bounds(q=16,P=40):
    B,D,tau=F(81),F(11757,10000),F(2,7)
    x=D*tau**(q-2); z=3*D*tau**(2*q-2)
    s=2*x/(1-tau); y=3*D*tau**(P-1)
    junction=B*(10+27/(1-tau))*(z+z*z)
    zoning=3*B*(9*s*(1+2*x)+2*(9*q+4)*(x+x*x))
    proxy=B*(27*(2*q+3)+1)*(y+y*y)
    return dict(junction=junction,zoning=zoning,proxy=proxy,total=junction+zoning)


def model_factory(q,P):
    R=BS.right_chains(max(P,2*q+3))
    common=left_chain(PE.end_vec(1,-1),P)[P]

    def incoming(t):
        return R[t[1]-1] if t[0]=='short' else R[P]

    @lru_cache(None)
    def arm(t,u,v):
        n=t[1] if t[0]=='short' else 2*q+t[1]
        J=mixed_junction(incoming(u),incoming(v))
        L=left_chain(J,n-1)
        if t[0]=='short':
            hw=[PE.hex_weights(L[i],R[n-1-i])[0] for i in range(n)]
            bw=[None]+[PE.base_weight(L[c-1],R[n-1-c]) for c in range(1,n)]
            leafL=L[-1]
        else:
            m=t[1]
            def pair(p,right_index):
                return (L[p] if p<q else common,
                        R[right_index] if p>=q+m else R[P])
            hw=[PE.hex_weights(*pair(i,n-1-i))[0] for i in range(n)]
            bw=[None]+[PE.base_weight(*pair(c-1,n-1-c)) for c in range(1,n)]
            leafL=common
        leaf=PE.dot(PE.vm(leafL,PE.W),PE.end_vec(1,1,1))/PE.dot(PE.vm(leafL,PE.W),PE.end_vec(1,1))
        data=(hw,bw,leaf)
        return data,BS.arm_states(*data)

    def model(types):
        h=[incoming(t) for t in types]
        J=mixed_junction(h[1],h[2]);Jminus=mixed_junction(h[1],h[2],power=2)
        central=PE.dot(PE.vm(Jminus,PE.W),h[0])/PE.dot(PE.vm(J,PE.W),h[0])
        arms=[];tables=[]
        for i,t in enumerate(types):
            u,v=sorted(types[j] for j in range(3) if j!=i)
            data,table=arm(t,u,v)
            arms.append(data);tables.append(table)
        best=None
        for states in product(*(sorted(table) for table in tables)):
            armcost=sum(table[s] for table,s in zip(tables,states))
            for C,T in BOOL:
                if not(C or T or any(s[1] for s in states)) or (C and all(s[2] for s in states)):
                    continue
                rankadd=3*int(C)+int(T and not C)+int(C and T and any(s[0]=='u' for s in states))
                deficit=armcost+rankadd-(3*int(C)+int(T))*central
                if best is None or deficit<best['deficit']:
                    best=dict(deficit=deficit,states=states,C=C,T=T)
        nsum=sum(t[1] if t[0]=='short' else 2*q+t[1] for t in types)
        fullcost=4*nsum-(sum(sum(w[2] for w in hw)+3*sum(bw[1:])+2*leaf for hw,bw,leaf in arms)+4*central)
        if all(t[0]=='short' for t in types):
            assert fullcost==0
        best.update(types=types,fullcost=fullcost,central=central)
        return best

    hw=PE.hex_weights(common,R[P])[0]
    bw=PE.base_weight(common,R[P])
    (gamma,_,_),trans=PA.gamma_cyc(PA.period_automaton(hw,bw))
    return model,gamma,trans


def uniform():
    OUT.mkdir(parents=True,exist_ok=True)
    q,P=16,40
    model,gamma,trans=model_factory(q,P)
    types=[('short',n) for n in range(2,2*q)]+[('long',m) for m in range(4)]
    records=[];worst={};counts={};maxfull=F(0)
    started=time.monotonic()
    for ts in combinations_with_replacement(types,3):
        row=model(ts)
        k=sum(t[0]=='long' for t in ts)
        counts[k]=counts.get(k,0)+1
        if k not in worst or row['deficit']<worst[k]['deficit']:
            worst[k]=row
        maxfull=max(maxfull,abs(row['fullcost']))
        records.append(row)
        if len(records)%250==0:
            print('models',len(records),'seconds',round(time.monotonic()-started,2),
                  'worst',float(min(r['deficit'] for r in worst.values())),flush=True)
    errs=error_bounds(q,P)
    modelmin=min(r['deficit'] for r in worst.values())
    full_bound=maxfull+errs['proxy']
    uniform_lower=min(modelmin-errs['proxy'],gamma-errs['proxy']-full_bound)-errs['total']
    summary=dict(q=q,proxy=P,counts=counts,worst=worst,model_minimum=modelmin,
                 max_abs_fullcost=maxfull,errors=errs,cycle_gap=gamma,full_bound=full_bound,
                 uniform_lower_bound=uniform_lower,
                 transitions={str(k):v[0] for k,v in trans.items()},
                 status='finite exact comparison; uniform implication requires written proof and review')
    (OUT/'uniform_cases.json').write_text(json.dumps(BS.serialize(records),indent=2)+'\n')
    (OUT/'uniform_summary.json').write_text(json.dumps(BS.serialize(summary),indent=2)+'\n')
    print('Counts',counts,'models',len(records),flush=True)
    for k,row in worst.items():
        print('Long arms',k,'minimum',float(row['deficit']),'at',row['types'],flush=True)
    print('Errors',{k:float(v) for k,v in errs.items()},flush=True)
    print('gamma',float(gamma),'max full',float(maxfull),'uniform bound',float(uniform_lower),flush=True)
    assert counts=={0:4960,1:1860,2:300,3:20}
    assert modelmin>F(29198,10**6)
    assert gamma>F(1259,10000)
    # The actual measured full cost already enters uniform_lower exactly.
    # This coarse display bound is sufficient for the printed proof.
    assert maxfull<F(1,10**15)
    assert uniform_lower>F(2698,100000)
    print('All declared exact comparisons PASS; written uniform argument remains necessary.',flush=True)


if __name__ == '__main__':
    started=time.monotonic()
    if '--uniform' in sys.argv:
        uniform()
    else:
        diagnostic()
    print('seconds',time.monotonic()-started,flush=True)
