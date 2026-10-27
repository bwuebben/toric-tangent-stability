#!/usr/bin/env python3
"""Exact finite verification for Local stability chambers for toric tangent bundles.

Requires Python 3.10+ and its standard library. No third-party modules.
Facet enumeration uses Cramer's rule; degrees use polar facet volumes;
flats use intersections of hyperplanes; matroid components use all circuits.
The mathematical justification of these operations is in the manuscript.
The exhaustive classification from which the 479 cases were selected is imported.
"""
from fractions import Fraction as F
from itertools import combinations
from functools import lru_cache
from collections import Counter
from pathlib import Path
import argparse
import json
import hashlib
import time

def det(a):
 a=[list(row) for row in a]; n=len(a); d=1; sign=1
 if not n:return 1
 for k in range(n-1):
  if a[k][k]==0:
   j=next((j for j in range(k+1,n) if a[j][k]),None)
   if j is None:return 0
   a[k],a[j]=a[j],a[k];sign=-sign
  pivot=a[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n): a[i][j]=(a[i][j]*pivot-a[i][k]*a[k][j])//d
  for i in range(k+1,n):a[i][k]=0
  d=pivot
 return sign*a[-1][-1]

def matrix_rank(a):
 a=[list(map(F,row)) for row in a]
 if not a:return 0
 row=0
 for col in range(len(a[0])):
  p=next((i for i in range(row,len(a)) if a[i][col]),None)
  if p is None:continue
  a[row],a[p]=a[p],a[row]
  for i in range(row+1,len(a)):
   if a[i][col]:
    q=a[i][col]/a[row][col]
    a[i]=[v-q*w for v,w in zip(a[i],a[row])]
  row+=1
  if row==len(a):break
 return row

def dot(a,b):return sum(x*y for x,y in zip(a,b))

def facet_enumeration(rays):
 n=len(rays[0]); facets={}; allids=range(len(rays))
 for ids in combinations(allids,n):
  A=[rays[i] for i in ids];d=det(A)
  if d==0:continue
  # Solve supporting equations exactly by Cramer's rule.
  normal=tuple(F(det([[1 if j==k else row[j] for j in range(n)] for row in A]),d) for k in range(n))
  if all(dot(u,normal)<=1 for u in rays):
   tight=tuple(i for i,u in enumerate(rays) if dot(u,normal)==1)
   assert len(tight)==n and abs(d)==1
   facets[tight]=normal
 return facets

def polar_facet_degrees(rays,facets):
 vertices=list(facets.values());m=len(rays)
 faces=[frozenset(j for j,v in enumerate(vertices) if dot(u,v)==1) for u in rays]
 @lru_cache(None)
 def triangulate(S):
  if len(S)==1:return [(next(iter(S)),)]
  apex=min(S)
  cuts={S & T for T in faces if S & T and S & T != S}
  boundary=[T for T in cuts if not any(T<U for U in cuts)]
  simplices=[]
  for T in boundary:
   if apex not in T:
    simplices += [(apex,)+t for t in triangulate(T)]
  return simplices
 result=[]
 for S in faces:
  simplexes=triangulate(S)
  assert all(len(t)==len(rays[0]) for t in simplexes)
  assert all(v.denominator==1 for x in vertices for v in x)
  result.append(sum(abs(det([[int(x) for x in vertices[j]] for j in t])) for t in simplexes))
 return result

def flats_via_hyperplanes(rays):
 n=len(rays[0]);m=len(rays);full=(1<<m)-1
 hyperplanes=set()
 for ids in combinations(range(m),n-1):
  A=[rays[i] for i in ids]
  normal=[(-1)**j*det([[row[k] for k in range(n) if k!=j] for row in A]) for j in range(n)]
  if not any(normal):continue
  mask=sum(1<<i for i,u in enumerate(rays) if dot(u,normal)==0)
  hyperplanes.add(mask)
 # A matroid flat is an intersection of its containing hyperplanes.
 flats={full}
 for h in hyperplanes:flats |= {f&h for f in list(flats)}
 flats.discard(full);flats.discard(0)
 return {f:matrix_rank([u for i,u in enumerate(rays) if f>>i&1]) for f in flats}

def components_by_circuits(rays):
 # Exhaust all minimal dependent ray subsets, independent of fundamental-circuit graph.
 m=len(rays);n=len(rays[0]);circuits=[]
 for k in range(2,n+2):
  for I in combinations(range(m),k):
   mask=sum(1<<i for i in I)
   if any(mask & c==c for c in circuits):continue
   if matrix_rank([rays[i] for i in I])<k:circuits.append(mask)
 groups=[1<<i for i in range(m)]
 for c in circuits:
  merged=c
  for g in groups:
   if g&c:merged |= g
  groups=[g for g in groups if not g&c]+[merged]
 return {frozenset(i for i in range(m) if g>>i&1) for g in groups}

def main():
    if not __debug__:
        raise SystemExit('Run without -O: mathematical checks require assertions.')
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root',type=Path,help='Optional complex_vb directory containing the original data/ and output/ files')
    args=parser.parse_args()
    folder=Path(__file__).resolve().parent
    data=json.loads((folder/'chambers.json').read_text())
    start=time.monotonic();counts=Counter();allflats=0
    for rec in data['records']:
     rays=rec['rays'];n=rec['dimension'];m=len(rays)
     facets=facet_enumeration(rays)
     assert set(facets)=={tuple(x) for x in rec['maximal_cones']},(rec['id'],'facets')
     degrees=polar_facet_degrees(rays,facets)
     assert degrees==rec['divisor_degrees'],(rec['id'],'degrees',degrees)
     flats=flats_via_hyperplanes(rays);allflats+=len(flats)
     assert len(flats)==rec['proper_flat_count'],(rec['id'],'flat count')
     total=sum(degrees);equal={}
     direction=list(map(F,rec['curve_class_direction']));limits=[]
     assert all(sum(rays[i][j]*direction[i] for i in range(m))==0 for j in range(n))
     for mask,rank in flats.items():
      margin=n*sum(d for i,d in enumerate(degrees) if mask>>i&1)-rank*total
      assert margin<=0
      derivative=n*sum(d for i,d in enumerate(direction) if mask>>i&1)-rank*sum(direction)
      if margin==0:
       equal[mask]=rank
       assert derivative<0
       assert all(sum(mask>>i&1 for i in f)==rank for f in facets)
      elif derivative>0:limits.append(F(-margin)/derivative)
     assert {sum(1<<i for i in e['rays']):e['rank'] for e in rec['equalities']}==equal
     expected_bound=min(limits) if limits else None
     assert expected_bound==(F(rec['strict_slope_interval_upper_bound']) if rec['strict_slope_interval_upper_bound'] else None)
     signatures={}
     for i in range(m):signatures.setdefault(tuple(bool(mask>>i&1) for mask in equal),[]).append(i)
     blocks=list(signatures.values())
     assert {tuple(b) for b in blocks}=={tuple(b) for b in rec['blocks']}
     B=rec['blocks'];s=len(B)
     C=[[int(i in b) for i in range(m)] for b in B]
     principal=list(map(list,zip(*rays)))
     assert matrix_rank(principal+C)-n==s,('C independence',rec['id'])
     ranks=[F(sum(degrees[i] for i in b)*n,total) for b in B]
     assert ranks==rec['block_ranks']
     assert len(equal)==2**(s-1)-1
     order={(p,q) for p in range(s) for q in range(s) if all(not(mask>>B[q][0]&1) or bool(mask>>B[p][0]&1) for mask in equal)}
     assert order=={(p,p) for p in range(s)}|{(0,q) for q in range(1,s)}
     counts[s-1]+=1
     print(rec['id'],'facets',len(facets),'flats',len(flats),'D^n',total,'degrees',degrees,'ranks',ranks,'elapsed',round(time.monotonic()-start,1),flush=True)
    assert len(data['records']) == 31 and allflats == 4890 and counts == {1:26,2:4,3:1}
    print('ALL 31 PASS; count',dict(counts),'proper flats',allflats,flush=True)

    products={r['id']:{frozenset(c) for c in r['components']} for r in data['product_decompositions']}
    indecomposable={r['id']:r for r in data['records']}
    selected={r['id']:r for r in json.loads((folder/'selected_rays.json').read_text())['records']}
    assert set(selected)==set(products)|set(indecomposable)
    assert len(selected)==479 and len(products)==448
    for n in range(3,7):
        subset={name:r for name,r in selected.items() if r['dimension']==n}
        for name,r in subset.items():
            rays=r['rays']
            if name in indecomposable:
                assert rays==indecomposable[name]['rays']
            expected=products.get(name,{frozenset(range(len(rays)))})
            assert components_by_circuits(rays)==expected,(name,'circuit components')
        assert len(subset)==data['selection'][str(n)]['strictly_semistable']
        print('DIMENSION',n,'all',len(subset),'selected circuit-component partitions PASS',flush=True)

    family=json.loads((folder/'projective_bundle_family.json').read_text())
    for identification in family['lattice_identifications']:
        rec=indecomposable[identification['id']]
        k=identification['k']; n=2*k
        basis=identification['basis_images']
        assert abs(det(basis))==1
        v=[0]*(n-1)+[1]
        rays=[v,[-x for x in v]]
        for i,a in enumerate([1]*k+[-1]*(k-1)):
            e=[int(j==i) for j in range(n)]
            w=[-x for x in e];w[-1]=a
            rays.extend([e,w])
        mapped=[tuple(sum(u[j]*basis[j][i] for j in range(n)) for i in range(n)) for u in rays]
        indices={tuple(u):i for i,u in enumerate(rec['rays'])}
        assert set(mapped)==set(indices)
        permutation=[indices[u] for u in mapped]
        assert permutation==identification['ray_indices_in_family_order']
        cones={tuple(sorted(permutation[2*i+((mask>>i)&1)] for i in range(n))) for mask in range(2**n)}
        assert cones=={tuple(c) for c in rec['maximal_cones']}
    print('Both projective-bundle fan identifications PASS',flush=True)

    if args.source_root is not None:
        base=args.source_root
        for source in data['sources']:
            assert hashlib.sha256((base/source['path']).read_bytes()).hexdigest()==source['sha256'],source['path']
        for n in range(3,7):
            raw=json.loads((base/f'data/smooth_toric_fano_{n}d.json').read_text())
            census=json.loads((base/f'output/sweep_{n}d.json').read_text())
            raw={r['id']:r['vertices'] for r in raw['polytopes']}
            ids={r['name'] for r in census if r['verdict']=='strictly semistable'}
            assert ids=={name for name,r in selected.items() if r['dimension']==n}
            assert all(raw[name]==selected[name]['rays'] for name in ids)
        print('Eight source hashes, census selections and all selected rays PASS',flush=True)
    else:
        print('Source-file comparison not requested; census selection remains imported.',flush=True)
    print('FINITE VERIFICATION COMPLETE; seconds',round(time.monotonic()-start,2),flush=True)

if __name__=='__main__':
    main()
