from decimal import Decimal as D,localcontext
import json
from pathlib import Path
with localcontext() as c:
 c.prec=70;h=D('1e-9');worst=0.;count=0
 for astr in ['0','0.21','0.3','0.49','0.7']:
  a=D(astr)
  for xs in ['0.1','0.4','0.9']:
   x=D(xs)
   def inverse(y):return 2*y/(1+a+((1+a)**2-4*a*y).sqrt())
   s=inverse(x);j=1+a-2*a*s
   for degree in [0,1,2]:
    def f(y):return inverse(y)**degree
    def power(k):
     if k>degree:return D(0)
     v=D(1)
     for q in range(k):v*=degree-q
     return v*s**(degree-k)
    exact=[power(1)/j,power(2)/j**2+power(1)*2*a/j**3,power(3)/j**3+3*power(2)*2*a/j**4+power(1)*12*a*a/j**5]
    fd=[(f(x+h)-f(x-h))/(2*h),(f(x+h)-2*f(x)+f(x-h))/h**2,(f(x+2*h)-2*f(x+h)+2*f(x-h)-f(x-2*h))/(2*h**3)]
    for l in range(3):
     err=float(abs(fd[l]-exact[l])/max(D(1),abs(exact[l])));worst=max(worst,err);count+=1
 assert worst<1e-12,worst
 out={'checks':count,'decimal_precision':70,'step':'1e-9','max_error_scaled_by_max_1_truth':worst,'method':'Independent centered finite differences of analytic inverse, orders 1, 2, 3; powers 0,1,2; parameters 0,.21,.3,.49,.7; x=.1,.4,.9'}
 p=Path(__file__).resolve().parents[1]/'measurements/derivative-accuracy/analytic_controls.json';p.write_text(json.dumps(out,indent=2));print(out)
