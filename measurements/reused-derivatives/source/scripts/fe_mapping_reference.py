"""Independent 90-digit explicit-inverse derivative references (no inverse jets)."""
from decimal import Decimal as D,localcontext
from math import comb,factorial
from functools import lru_cache
from itertools import product

def inverse_power_derivatives(x,a,power,order):
 if not a:
  return [D(factorial(power)//factorial(power-k))*x**(power-k) if k<=power else D(0) for k in range(order+1)]
 disc=(1+a)**2-4*a*x;root=disc.sqrt();s=[2*x/(1+a+root)]
 for k in range(1,order+1):
  odd=1
  for j in range(1,2*k-2,2):odd*=j
  s.append(D(odd)*(2*a)**(k-1)/root**(2*k-1))
 out=[D(1)]+[D(0)]*order
 for _ in range(power):out=[sum(D(comb(k,j))*out[j]*s[k-j] for j in range(k+1)) for k in range(order+1)]
 return out

def add(*polys):
 out={}
 for p in polys:
  for k,v in p.items():out[k]=out.get(k,D(0))+v
 return {k:v for k,v in out.items() if v}
def scale(a,p):return {k:a*v for k,v in p.items()}
def mul(a,b):
 out={}
 for i,v in a.items():
  for j,w in b.items():
   k=tuple(x+y for x,y in zip(i,j));out[k]=out.get(k,D(0))+v*w
 return {k:v for k,v in out.items() if v}
def power(a,n):
 v={(0,0,0):D(1)}
 for _ in range(n):v=mul(v,a)
 return v
@lru_cache(None)
def shear_polynomial(curved):
 x={(1,0,0):D(1)};y={(0,1,0):D(1)};z={(0,0,1):D(1)}
 t=add(y,scale(D('-.15'),z),scale(D('.15'),mul(z,z))) if curved else y
 s=add(x,scale(D('-.2'),t),scale(D('.2'),mul(t,t))) if curved else x
 return add(power(s,5),power(t,5),power(z,5),mul(power(s,2),power(t,2)),mul(mul(s,t),z))
@lru_cache(None)
def derivatives_at(dim,curved,x,y,z,precision=90):
 with localcontext() as ctx:
  ctx.prec=precision;x,y,z=map(D,(x,y,z));out={}
  if dim==2:
   sx={p:inverse_power_derivatives(x,D('.7') if curved else D(0),p,10) for p in [0,2,5]}
   ty={p:inverse_power_derivatives(y,D('.49') if curved else D(0),p,10) for p in [0,2,5]}
   for nx in range(11):
    for ny in range(11-nx):out[nx,ny,0]=sx[5][nx]*ty[0][ny]+sx[0][nx]*ty[5][ny]+sx[2][nx]*ty[2][ny]
  else:
   poly=shear_polynomial(curved);coords=(x,y,z)
   for nx in range(11):
    for ny in range(11-nx):
     for nz in range(11-nx-ny):
      alpha=(nx,ny,nz);v=D(0)
      for p,c in poly.items():
       if any(p[d]<alpha[d] for d in range(3)):continue
       term=c
       for d in range(3):term*=D(factorial(p[d])//factorial(p[d]-alpha[d]))*coords[d]**(p[d]-alpha[d])
       v+=term
      out[alpha]=v
  return out

def explicit_value(dim,curved,coords):
 x,y=coords[:2]
 if dim==2:
  s=inverse_power_derivatives(x,D('.7') if curved else D(0),1,0)[0]
  t=inverse_power_derivatives(y,D('.49') if curved else D(0),1,0)[0]
  return s**5+t**5+s*s*t*t
 z=coords[2];t=y-D('.15')*z*(1-z) if curved else y;s=x-D('.2')*t*(1-t) if curved else x
 return s**5+t**5+z**5+s*s*t*t+s*t*z

def finite_difference(dim,curved,coords,alpha,step):
 with localcontext() as ctx:
  ctx.prec=90;h=D(step);coords=list(map(D,coords));stencils={0:[(0,D(1))],1:[(-1,D('-.5')),(1,D('.5'))],2:[(-1,D(1)),(0,D(-2)),(1,D(1))],3:[(-2,D('-.5')),(-1,D(1)),(1,D(-1)),(2,D('.5'))]}
  val=D(0)
  for terms in product(*(stencils[k] for k in alpha[:dim])):
   p=[coords[d]+terms[d][0]*h for d in range(dim)];w=D(1)
   for _,a in terms:w*=a
   val+=w*explicit_value(dim,curved,p)
  return val/h**sum(alpha)
