"""Fixed E1 end-attachment experiment; new files only, exact final evidence.

Floating-point column generation proposes bases; all final decompositions,
exchange graphs, degrees and negative witnesses are checked over Q.
The candidate list is declared in research_notes/Toric interactions execution/
end_attachments.md. No search beyond that list occurs here.
"""
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import comb, factorial
from pathlib import Path
import json
import sys
import time
sys.set_int_max_str_digits(100000)
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / 'src'))
# Portable theorem computation: numerical basis search is not included.
from types import SimpleNamespace
from toric_stability import rref, in_span, rank as exact_rank
MC = SimpleNamespace(rref=rref, in_span=in_span, exact_rank=exact_rank)
import verify_certificate as VC
import tree_family as TF
import prism_check as PC

OUT = HERE.parent / 'output/research_end_attachments'
HEX = [(1,0),(0,1),(1,1),(-1,0),(0,-1),(-1,-1)]
QVERT = [(1,0),(1,-1),(0,-1),(-1,0),(-1,1),(0,1)]
CAPS = {0:[(-1,0)], 1:[(0,1),(-1,-1)],
        2:[(-1,0),(0,1),(0,-1)], 3:[(1,0),(-1,0),(0,1),(-1,-1)]}
PATTERNS = [(1,0),(2,0),(1,1),(3,0),(1,2)]

@lru_cache(None)
def moment(a,b,facet=None):
    if facet is not None:
        vs=[v for v in QVERT if sum(x*y for x,y in zip(facet,v)) == -1]
        assert len(vs)==2
        (x,y),(xx,yy)=vs; dx,dy=xx-x,yy-y
        return sum(F(comb(a,i)*comb(b,k)*x**(a-i)*y**(b-k)*dx**i*dy**k,i+k+1)
                   for i in range(a+1) for k in range(b+1))
    ans=F(0)
    for u,v in zip(QVERT,QVERT[1:]+QVERT[:1]):
        det=abs(u[0]*v[1]-u[1]*v[0])
        for i in range(a+1):
            for k in range(b+1):
                ans += F(det*comb(a,i)*comb(b,k)*u[0]**i*v[0]**(a-i)*u[1]**k*v[1]**(b-k)
                         *factorial(i+k)*factorial(a+b-i-k),factorial(a+b+2))
    return ans

class Attachment:
    def __init__(self,j,left,right):
        self.j,self.left,self.right=j,left,right
        self.name=f'A{left}_A{right}_j{j}'
        self.tw=[]; self.b=[]
        for t in CAPS[left]:
            self.tw.append({0:t}); self.b.append(1)
        for i in range(j-1):
            self.tw.append({i:(1,0),i+1:(-1,0)}); self.b.append(2)
        for t in CAPS[right]:
            self.tw.append({j-1:tuple(-x for x in t)}); self.b.append(1)
        self.k=len(self.b); self.n=2*j+sum(self.b)
        self.rays=[]; self.labels=[]; off=2*j
        for i in range(j):
            for h in HEX:
                v=[0]*self.n;v[2*i:2*i+2]=h
                self.rays.append(tuple(v));self.labels.append(('hex',i,h))
        for c,(b,tw) in enumerate(zip(self.b,self.tw)):
            for a in range(b):
                v=[0]*self.n;v[off+a]=1
                self.rays.append(tuple(v));self.labels.append(('base',c,a+1))
            v=[0]*self.n
            for a in range(b):v[off+a]=-1
            for i,t in tw.items():v[2*i:2*i+2]=t
            self.rays.append(tuple(v));self.labels.append(('base',c,0))
            off+=b
        assert self.n==4*j+left+right
        assert len(self.rays)-self.n==5*j+left+right+1
        assert all(len(tw)==b and all(t in HEX for t in tw.values()) for b,tw in zip(self.b,self.tw))

    def poly(self,drop=None):
        ans={(0,)*(2*self.j):F(1)}
        for c,(b,tw) in enumerate(zip(self.b,self.tw)):
            power=b-int(c==drop)
            for _ in range(power):
                nxt={}
                for mo,coef in ans.items():
                    nxt[mo]=nxt.get(mo,0)+(b+1)*coef
                    for i,t in tw.items():
                        for z,val in enumerate(t):
                            if val:
                                mm=list(mo);mm[2*i+z]+=1;mm=tuple(mm)
                                nxt[mm]=nxt.get(mm,0)+val*coef
                ans={mo:coef for mo,coef in nxt.items() if coef}
            if power>1:ans={mo:coef/F(factorial(power)) for mo,coef in ans.items()}
        return ans

    def integrate(self,poly,facet=None):
        ans=F(0)
        for mo,coef in poly.items():
            for i in range(self.j):
                coef*=moment(mo[2*i],mo[2*i+1],facet[1] if facet and facet[0]==i else None)
                if not coef:break
            ans+=coef
        return ans

    def degrees(self):
        poly=self.poly(); volume=self.integrate(poly)
        hd={(i,h):self.integrate(poly,(i,h))/volume for i in range(self.j) for h in HEX}
        bd=[self.integrate(self.poly(c))/volume for c in range(self.k)]
        delta=[hd[(l[1],l[2])] if l[0]=='hex' else bd[l[1]] for l in self.labels]
        assert sum(delta)==self.n and all(x>0 for x in delta)
        # Principal divisor degree balance in every lattice direction.
        assert all(sum(d*v[a] for d,v in zip(delta,self.rays))==0 for a in range(self.n))
        return delta,volume*factorial(self.n)


def controls():
    for a in range(7):
        for b in range(7):
            assert moment(a,b)==PC._int_hex_mono(a,b)
            for h in HEX:assert moment(a,b,h)==PC.int_edge({(a,b):F(1)},h)
    for j in [1,2]:
        obj=Attachment(j,0,0);delta,D=obj.degrees();old=TF.path(j);deg,D0=old.degrees()
        assert obj.rays==old.rays and D==D0 and delta==[old.n*d/D0 for d in deg]
    for left,right in PATTERNS:
        obj=Attachment(1,left,right);delta,D=obj.degrees()
        twists=[tw[0] for tw in obj.tw]
        assert all(sum(t[a] for t in twists)==0 for a in range(2))
        hd,bd,D0,mu,poly=PC.prism_degrees(twists)
        expected=[hd[h]/mu for h in HEX]+[d/mu for d in bd for _ in range(2)]
        assert D==D0 and delta==expected
    print('Controls passed: 343 independent polygon/facet moments, ordinary paths j=1,2, and all five zero-sum j=1 degree controls.',flush=True)


def run(j,left,right):
    obj=Attachment(j,left,right);start=time.time()
    delta,D=obj.degrees()
    print(f'{obj.name}: n={obj.n}, rho={len(obj.rays)-obj.n}; exact degrees done ({time.time()-start:.2f}s)',flush=True)
    cert=MC.certify(obj.rays,delta,name=obj.name,tries=4,log=lambda _:None)
    rec=dict(name=obj.name,j=j,left=left,right=right,n=obj.n,rho=len(obj.rays)-obj.n,
             base_dimensions=obj.b,twists=[{str(i):list(t) for i,t in tw.items()} for tw in obj.tw],
             rays=obj.rays,labels=obj.labels,delta=list(map(str,delta)),anticanonical_degree=str(D))
    if cert['status'].startswith('stable'):
        path=OUT/(obj.name+'_certificate.json')
        MC.save(cert,obj.rays,delta,path,'research_end_attachments.py: exact polynomial simplex slices; independently checked polygon moments and j=1 prism degrees')
        VC.verify(path)
        rec.update(status='stable',certificate=path.name,bases=len(cert['bases']))
        print(f'  STABLE: {len(cert["bases"])} exact rational bases; standalone verifier passed.',flush=True)
    else:
        violation=cert.get('violation')
        if violation is not None:
            amount,rank,indices=violation
            assert amount==sum(delta[q] for q in indices)-rank and amount>=0
            rk,canon=MC.rref([obj.rays[q] for q in indices]);assert rk==rank
            closure=[q for q,v in enumerate(obj.rays) if MC.in_span(canon,v)];assert closure==indices
            rec.update(status='unstable' if amount>0 else 'not stable (equality)',witness=dict(
                rank=rank,ray_indices=indices,labels=[obj.labels[q] for q in indices],
                excess_weight=str(amount),basis=[list(map(str,row)) for row in canon]))
            print(f'  {rec["status"].upper()}: rank {rank}, excess normalized weight {amount}; {rec["witness"]["labels"]}',flush=True)
        else:
            rec['status']='unresolved: no exact certificate or witness'
            print('  UNRESOLVED: numerical search supplied no exact conclusion.',flush=True)
    rec['seconds']=time.time()-start
    (OUT/(obj.name+'.json')).write_text(json.dumps(rec,indent=2)+'\n')
    return rec

def finite_main():
    OUT.mkdir(parents=True,exist_ok=True)
    controls()
    results=[]
    for j in [1,2,3,4]:
        for left,right in PATTERNS:
            results.append(run(j,left,right))
            (OUT/'summary.json').write_text(json.dumps([{k:r[k] for k in ['name','j','left','right','n','rho','status','seconds']} for r in results],indent=2)+'\n')

# Uniform proof: finite local end pieces and the established path kernel.
import path_transfer as PT
import path_all_j as PA
import path_flats as PF

def pmul2(P,Q):
    out={}
    for a,x in P.items():
        for b,y in Q.items():
            k=(a[0]+b[0],a[1]+b[1]);out[k]=out.get(k,0)+x*y
    return {k:v for k,v in out.items() if v}

@lru_cache(None)
def cap_vectors(a):
    def polynomial(drop=None):
        out={(0,0):F(1)}
        for i,t in enumerate(CAPS[a]):
            if i==drop:continue
            fac={(0,0):F(2)}
            if t[0]:fac[(1,0)]=F(t[0])
            if t[1]:fac[(0,1)]=F(t[1])
            out=pmul2(out,fac)
        return out
    def vec(P,h=None):
        return tuple(sum(v*moment(x+p,y,h) for (x,y),v in P.items()) for p in range(3))
    poly=polynomial(); bulk=vec(poly)
    modified=[vec(poly,h) for h in HEX]
    for i in range(len(CAPS[a])):modified.extend([vec(polynomial(i))]*2)
    return bulk,tuple(modified)

def cap_rays(a):
    n=2+len(CAPS[a]); rays=[]
    for h in HEX:rays.append(tuple(h)+(0,)*(n-2))
    for i,t in enumerate(CAPS[a]):
        v=[0]*n;v[2+i]=1;rays.append(tuple(v))
        v=[0]*n;v[:2]=t;v[2+i]=-1;rays.append(tuple(v))
    return rays

@lru_cache(None)
def cap_flats(a):
    rays=cap_rays(a); n=len(rays[0]); root=rays[0]
    states={():dict(indices=(),rank=0,ground=False,quotient_full=False)}; todo=[()]
    while todo:
        idx=todo.pop(); base=[rays[q] for q in idx]
        for e in range(len(rays)):
            if e in idx:continue
            rank,canon=MC.rref(base+[rays[e]])
            closed=tuple(q for q,v in enumerate(rays) if MC.in_span(canon,v))
            if closed in states:continue
            ground=MC.in_span(canon,root)
            states[closed]=dict(indices=closed,rank=rank,ground=ground,
                               quotient_full=(rank+int(not ground)==n))
            todo.append(closed)
    result=tuple(states.values())
    print(f'Cap A{a}: {len(result)} exact local ray flats in rank {n}.',flush=True)
    return result


def transfer_vectors(j,left,right):
    mL,_=cap_vectors(left); mR0,_=cap_vectors(right)
    mR=tuple((-1)**p*x for p,x in enumerate(mR0))
    L=[list(mL)]; R=[list(mR)]
    for _ in range(j-2):
        L.append(PT.vm(PT.vm(L[-1],PT.KAP),PT.W))
        R.append(PT.mv(PT.W,PT.mv(PT.KAP,R[-1])))
    return L,R


def transfer_weights(j,left,right):
    assert j>=2
    L,R=transfer_vectors(j,left,right)
    den=PT.dot(PT.vm(L[0],PT.KAP),R[j-2])
    _,modsL=cap_vectors(left);_,modsR=cap_vectors(right)
    capL=[PT.dot(PT.vm(v,PT.KAP),R[j-2])/den for v in modsL]
    capR=[PT.dot(PT.vm(L[j-2],PT.KAP),[(-1)**p*x for p,x in enumerate(v)])/den for v in modsR]
    # Mirror identifies the right cap's local h ray with the global -h ray.
    hw={}
    for i in range(1,j-1):
        hw[i]=PA.hex_weights(PT.vm(L[i-1],PT.KAP),PT.mv(PT.KAP,R[j-2-i]))
    bw={c:PT.dot(PT.vm(L[c-1],PT.KAP1),R[j-1-c])/den for c in range(1,j)}
    return capL,capR,hw,bw,den


def local_step(st,C,Phi,h):
    link=Phi and C
    rank=2*int(C)+int(Phi and not C)
    if h:
        rank+=h+int(link and st=='u'); ns='G'
    elif link:rank+=1;ns='g' if st!='u' else 'u'
    else:ns='u'
    return rank,ns


def end_minimum(j,left,right,weights):
    capL,capR,hw,bw,_=weights
    # qfull: adjoining only the outgoing distinguished line gives full rank.
    cur={}
    def add(out,key,value,hist):
        if key not in out or value<out[key][0]:out[key]=(value,hist)
    for opt in cap_flats(left):
        cost=F(opt['rank'])-sum(capL[q] for q in opt['indices'])
        key=('G' if opt['ground'] else 'u',bool(opt['indices']),opt['quotient_full'])
        add(cur,key,cost,(opt['indices'],))
    for i in range(1,j-1):
        nxt={}
        for (st,nonempty,qfull),(cost,hist) in cur.items():
            for C,Phi in PF.BOOL2:
                for h in PF.HEXOPT:
                    inc,ns=local_step(st,C,Phi,h)
                    full=qfull and (int(st=='u')+4-inc==int(ns=='u'))
                    value=cost+inc-bw[i]*(2*C+Phi)-hw[i][h]
                    add(nxt,(ns,nonempty or C or Phi or h>0,full),value,hist+((C,Phi,h),))
        cur=nxt
    best=None;dimR=len(cap_rays(right)[0])
    for (st,nonempty,qfull),(cost,hist) in cur.items():
        for C,Phi in PF.BOOL2:
            for opt in cap_flats(right):
                inc=2*int(C)+int(Phi and not C)+opt['rank']
                if C and Phi:inc+=int(st=='u' or not opt['ground'])
                full=qfull and int(st=='u')+2+dimR-inc==0
                if full or not(nonempty or C or Phi or opt['indices']):continue
                value=cost+inc-bw[j-1]*(2*C+Phi)-sum(capR[q] for q in opt['indices'])
                if best is None or value<best[0]:best=(value,hist+((C,Phi),opt['indices']))
    return best


def transfer_gate():
    for j in [2,3,4]:
        for left,right in PATTERNS:
            obj=Attachment(j,left,right);d,D=obj.degrees()
            capL,capR,hw,bw,Z=transfer_weights(j,left,right)
            assert factorial(obj.n)*Z==D
            # Exact raywise alignment, including the reflected right cap.
            left_count=len(CAPS[left]);first_right=left_count+j-1
            expected=[]
            for lab in obj.labels:
                if lab[0]=='hex':
                    i,h=lab[1:]
                    if i==0:value=capL[HEX.index(h)]
                    elif i==j-1:value=capR[HEX.index(tuple(-x for x in h))]
                    else:
                        L,R=transfer_vectors(j,left,right)
                        lv,rv=PT.vm(L[i-1],PT.KAP),PT.mv(PT.KAP,R[j-2-i])
                        value=PT.dot(PT.vm(lv,PT.WH[PT.KIND[h]]),rv)/PT.dot(PT.vm(lv,PT.W),rv)
                else:
                    c=lab[1]
                    if c<left_count:value=capL[6+2*c]
                    elif c>=first_right:value=capR[6+2*(c-first_right)]
                    else:value=bw[c-left_count+1]
                expected.append(value)
            assert d==expected,(j,left,right)
    print('Transfer formulas agree with every polynomial-derived ray degree at j=2,3,4, all five patterns.',flush=True)


def config_indices(obj,hist):
    left=hist[0];right=hist[-1];j=obj.j;L=len(CAPS[obj.left]);Rstart=L+j-1
    index={lab:q for q,lab in enumerate(obj.labels)}
    out=[]
    for q in left:
        if q<6:out.append(index[('hex',0,HEX[q])])
        else:
            leaf=(q-6)//2;typ=1 if(q-6)%2==0 else 0
            out.append(index[('base',leaf,typ)])
    for i,opt in enumerate(hist[1:-2],start=1):
        C,Phi,h=opt
        if C:out.extend(index[('base',L+i-1,a)] for a in [1,2])
        if Phi:out.append(index[('base',L+i-1,0)])
        if h:out.extend(index[('hex',i,hh)] for hh in ([(1,0),(-1,0)] if h==1 else HEX))
    C,Phi=hist[-2]
    if C:out.extend(index[('base',Rstart-1,a)] for a in [1,2])
    if Phi:out.append(index[('base',Rstart-1,0)])
    for q in right:
        if q<6:out.append(index[('hex',j-1,tuple(-x for x in HEX[q]))])
        else:
            leaf=(q-6)//2;typ=1 if(q-6)%2==0 else 0
            out.append(index[('base',Rstart+leaf,typ)])
    return sorted(out)


def cap_model_weights(m,q,KP,left,right):
    j=2*q+m
    # Arbitrary endpoints converge to the same bulk functions after propagation.
    L,R=transfer_vectors(KP+3,left,right)
    capL=[PT.dot(PT.vm(v,PT.KAP),R[KP])/PT.dot(PT.vm(L[0],PT.KAP),R[KP]) for v in cap_vectors(left)[1]]
    rm=[(-1)**p*x for p,x in enumerate(cap_vectors(right)[0])]
    den=PT.dot(PT.vm(L[KP],PT.KAP),rm)
    capR=[PT.dot(PT.vm(L[KP],PT.KAP),[(-1)**p*x for p,x in enumerate(v)])/den for v in cap_vectors(right)[1]]
    hw={};bw={}
    def zone(pos):return 'left' if pos<q else ('deep' if pos<q+m else 'right')
    for i in range(1,j-1):
        z=zone(i)
        lv=PT.vm(L[i-1] if z=='left' else L[KP],PT.KAP)
        rv=PT.mv(PT.KAP,R[j-2-i] if z=='right' else R[KP])
        hw[i]=PA.hex_weights(lv,rv)
    for c in range(1,j):
        z=zone(c)
        lv=L[c-1] if z=='left' else L[KP]
        rv=R[j-1-c] if z=='right' else R[KP]
        bw[c]=PT.dot(PT.vm(lv,PT.KAP1),rv)/PT.dot(PT.vm(lv,PT.KAP),rv)
    return capL,capR,hw,bw,None


def uniform_main():
    OUT.mkdir(parents=True,exist_ok=True)
    transfer_gate()
    for a in range(4):
        flats=cap_flats(a)
        (OUT/f'cap_A{a}_flats.json').write_text(json.dumps(list(flats),indent=2)+'\n')
    import random
    rng=random.Random(1709)
    for j in [2,3,4]:
        for left,right in PATTERNS:
            obj=Attachment(j,left,right);d,D=obj.degrees();weights=transfer_weights(j,left,right)
            val,hist=end_minimum(j,left,right,weights)
            idx=config_indices(obj,hist);rk=MC.exact_rank([obj.rays[q] for q in idx])
            assert 0<rk<obj.n and val==rk-sum(d[q] for q in idx)
            # Random full configurations test rank and full-space flags independently.
            for _ in range(200):
                lo=rng.choice(cap_flats(left));ro=rng.choice(cap_flats(right))
                hist=(lo['indices'],)+tuple((rng.choice([False,True]),rng.choice([False,True]),rng.randrange(3)) for _ in range(j-2))
                hist+=( (rng.choice([False,True]),rng.choice([False,True])),ro['indices'])
                st='G' if lo['ground'] else 'u'; rank=lo['rank'];qfull=lo['quotient_full']
                for C,Phi,h in hist[1:-2]:
                    inc,ns=local_step(st,C,Phi,h)
                    qfull=qfull and int(st=='u')+4-inc==int(ns=='u');rank+=inc;st=ns
                C,Phi=hist[-2]
                inc=2*int(C)+int(Phi and not C)+ro['rank']+int(C and Phi and(st=='u' or not ro['ground']))
                full=qfull and int(st=='u')+2+len(cap_rays(right)[0])-inc==0
                rank+=inc
                idx=config_indices(obj,hist);exact=MC.exact_rank([obj.rays[q] for q in idx]) if idx else 0
                assert rank==exact and full==(exact==obj.n),(j,left,right,hist,rank,exact,full)
            print(f'End-state ranks checked: {obj.name}, minimum {float(val):.12g}',flush=True)
    q,KP=20,40
    results=[]
    for left,right in PATTERNS:
        print(f'Uniform pattern A{left}/A{right}:',flush=True)
        minima=[];fullerrors=[]
        for m in range(4):
            weights=cap_model_weights(m,q,KP,left,right)
            val,hist=end_minimum(2*q+m,left,right,weights)
            capL,capR,hw,bw,_=weights
            n=4*(2*q+m)+left+right
            fullerror=abs(F(n)-sum(capL)-sum(capR)-sum(x[2] for x in hw.values())-3*sum(bw.values()))
            minima.append(val);fullerrors.append(fullerror)
            print(f'  m={m}: model minimum {float(val):.12g}, full error {float(fullerror):.3e}',flush=True)
        finite=[]
        for j in range(2,2*q+5):
            val,hist=end_minimum(j,left,right,transfer_weights(j,left,right))
            assert val>0,(left,right,j,val)
            finite.append((val,j))
        bulkL,bulkR=PA.chain_vectors(KP)
        hw=PA.hex_weights(bulkL[KP],bulkR[KP]);bw=PA.base_weight(bulkL[KP],bulkR[KP])
        gamma=PA.gamma_cyc(PA.period_automaton(hw,bw))[0][0]
        tau,Delta,B=F(2,7),F(11757,10000),F(243)
        x=Delta*tau**(q-3);y=2*Delta*tau**(KP-2)
        deep=B*9*(2*x/(1-tau))*(1+2*x)
        zones=B*2*(9*q+10)*(x+x*x)
        proxy=B*(9*(2*q+4)+10)*(y+y*y)
        budget=deep+zones+2*proxy+max(fullerrors)
        bound=min([gamma]+minima)-budget
        assert bound>0
        print(f'  finite j=2..{2*q+4} positive, smallest {float(min(finite)[0]):.12g} at j={min(finite)[1]}; budget {float(budget):.3e}; lower bound {float(bound):.12g}',flush=True)
        results.append(dict(left=left,right=right,q=q,proxy_index=KP,
            model_minima=list(map(str,minima)),full_errors=list(map(str,fullerrors)),
            finite_minimum=str(min(finite)[0]),finite_minimum_j=min(finite)[1],
            gamma=str(gamma),budget=str(budget),uniform_lower_bound=str(bound)))
        (OUT/'uniform_proposal.json').write_text(json.dumps(results,indent=2)+'\n')

if __name__=='__main__':
    if '--uniform' in sys.argv:uniform_main()
    else:raise SystemExit('Run with --uniform; optional numerical basis search remains in the research source.')
