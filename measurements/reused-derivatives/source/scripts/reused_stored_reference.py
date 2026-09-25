"""Independent Decimal interpolation, Newton inversion and physical differences.
Only stored binary64 nodes/coefficients enter the reference. No chain-rule
recurrence or inverse Taylor series from either benchmark is reused.
"""
from decimal import Decimal as D,localcontext
from fractions import Fraction as F
from functools import lru_cache
from itertools import product
from math import factorial
import csv,json,time
from fe_mapping_reference import derivatives_at
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'measurements/reused-derivatives'
def solve(a,b):
 a=[list(row)+[v] for row,v in zip(a,b)];n=len(b)
 for i in range(n):
  pivot=max(range(i,n),key=lambda k:abs(a[k][i]));a[i],a[pivot]=a[pivot],a[i];v=a[i][i];a[i]=[x/v for x in a[i]]
  for j in range(n):
   if i==j:continue
   v=a[j][i];a[j]=[x-v*y for x,y in zip(a[j],a[i])]
 return [a[i][-1] for i in range(n)]
def mul(a,b):
 c=[D(0)]*(len(a)+len(b)-1)
 for i,x in enumerate(a):
  for j,y in enumerate(b):c[i+j]+=x*y
 return c
class Polynomial:
 def __init__(self,rows,dim,component):
  axes=[sorted({r[d] for r in rows}) for d in range(dim)];polys=[]
  for axis in axes:
   ps=[]
   for x in axis:
    p=[D(1)];den=D(1)
    for y in axis:
     if x!=y:p=mul(p,[-y,D(1)]);den*=x-y
    ps.append([v/den for v in p])
   polys.append(ps)
  keys=list(product(*(range(len(a)) for a in axes)));table={tuple(axes[d].index(row[d]) for d in range(dim)):row[3+component] for row in rows}
  for d in range(dim):
   nxt={}
   for k in keys:nxt[k]=sum(polys[d][j][k[d]]*table[k[:d]+(j,)+k[d+1:]] for j in range(len(axes[d])))
   table=nxt
  self.terms=table;self.dim=dim;self.degree=[len(a)-1 for a in axes]
 def value(self,x,derivative=-1):
  powers=[[D(1)] for _ in x]
  for d in range(self.dim):
   for k in range(self.degree[d]):powers[d].append(powers[d][-1]*x[d])
  out=D(0)
  for k,v in self.terms.items():
   if derivative>=0:
    if not k[derivative]:continue
    v*=k[derivative]
   for d in range(self.dim):v*=powers[d][k[d]-(d==derivative)]
   out+=v
  return out
@lru_cache(None)
def stencil(order):
 if not order:return [(0,D(1))]
 radius=(order+4)//2;xs=list(range(-radius,radius+1));n=len(xs)
 weights=solve([[F(x)**k for x in xs] for k in range(n)],[F(factorial(order)) if k==order else F(0) for k in range(n)])
 return [(x,D(w.numerator)/D(w.denominator)) for x,w in zip(xs,weights) if w]
def fixture(path):
 lines=iter(path.read_text().splitlines());dim,n,curved=map(int,next(lines).split());fields=[]
 for which in range(2):
  count,components=map(int,next(lines).split());rows=[[D.from_float(float(x)) for x in next(lines).split()] for _ in range(count)]
  fields.append([Polynomial(rows,dim,c) for c in range(components)])
 return dim,n,fields[0],fields[1][1]
def main():
 start=time.monotonic();raw=list(csv.DictReader((OUT/'high-order.csv').open()));results=[]
 with localcontext() as ctx:
  ctx.prec=90
  for dim in [2,3]:
   for n in [8,32]:
    _,_,geo,u=fixture(OUT/'fixtures'/f'high-{dim}-{n}-1-10.txt');xi=[D.from_float(float(v)) for v in ['.23','.31','.37'][:dim]];x=[g.value(xi) for g in geo]
    def evaluate(physical):
     r=xi[:]
     for iteration in range(12):
      residual=[geo[d].value(r)-physical[d] for d in range(dim)]
      if max(map(abs,residual))<D('1e-80'):break
      delta=solve([[g.value(r,d) for d in range(dim)] for g in geo],residual);r=[v-w for v,w in zip(r,delta)]
     else:raise RuntimeError('Newton did not converge')
     return u.value(r)
    alphas=[(1,0,0),(3,0,0),(5,0,0),(2,1,0),(3,2,0)] if dim==2 else [(1,0,0),(0,0,3),(0,0,5),(1,1,1),(1,1,3)]
    for alpha in alphas:
     refs=[]
     for step in ['.0001','.00005']:
      h=D(step)/n;value=D(0)
      for terms in product(*(stencil(k) for k in alpha[:dim])):
       coords=[x[d]+h*terms[d][0] for d in range(dim)];weight=D(1)
       for _,w in terms:weight*=w
       value+=weight*evaluate(coords)
      refs.append(value/h**sum(alpha))
     sensitivity=abs(refs[1]-refs[0])/max(D(1),abs(refs[1]));assert sensitivity<D('1e-7'),(dim,n,alpha,sensitivity)
     matches=[r for r in raw if (int(r['dimension']),int(r['n']),r['curved'],r['capacity'],r['point'],r['field'])==(dim,n,'1','10','0','quintic') and tuple(int(r[k]) for k in ['dx','dy','dz'])==alpha]
     ideal=derivatives_at(dim,1,*([str(v) for v in x]+(['0'] if dim==2 else [])))[alpha]
     for r in matches:
      value=D.from_float(float(r['value']));error=abs(value-refs[1]);results.append(dict(dimension=dim,n=n,dx=alpha[0],dy=alpha[1],dz=alpha[2],order=sum(alpha),method=r['method'],value=str(value),stored_reference=str(refs[1]),ideal_reference=str(ideal),representation_error=str(abs(refs[1]-ideal)),total_error=str(abs(value-ideal)),element_scaled_error=str(error/D(n)**sum(alpha)),absolute_error=str(error),scaled_error=str(error/max(D(1),abs(refs[1]))),relative_error=str(error/abs(refs[1])) if refs[1] else 'nan',step_sensitivity=str(sensitivity)))
     print('stored reference',dim,n,alpha,'step sensitivity',float(sensitivity),flush=True)
 with (OUT/'stored-reference.csv').open('w') as f:w=csv.DictWriter(f,fieldnames=results[0]);w.writeheader();w.writerows(results)
 (OUT/'stored-reference-run.json').write_text(json.dumps(dict(seconds=time.monotonic()-start,precision=90,steps=['h*1e-4','h*5e-5'],probes=len(results)//5),indent=2))
if __name__=='__main__':main()
