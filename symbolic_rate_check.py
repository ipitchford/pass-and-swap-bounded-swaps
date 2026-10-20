"""Verify the two product-form obstructions as identities in one indeterminate.

Only Python integers and fractions are used. Polynomial equality is decided by
cross multiplication, so no sampled parameter value is used as a proof.
"""
from fractions import Fraction as F
from itertools import permutations
from pathlib import Path
import json
from independent_check import balanced, orientation, transitions

def check(value,message):
    if not value:raise RuntimeError(message)

def trim(a):
    a=list(a)
    while len(a)>1 and a[-1]==0:a.pop()
    return tuple(a)

def add(a,b):
    n=max(len(a),len(b))
    return trim([(a[i] if i<len(a) else 0)+(b[i] if i<len(b) else 0) for i in range(n)])

def mul(a,b):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    return trim(c)

class RF:
    def __init__(self,n=(0,),d=(1,)):
        if isinstance(n,int):n=(n,)
        self.n,self.d=trim(n),trim(d)
        check(self.d!=(0,),'zero polynomial denominator')
    def __add__(self,x):
        x=as_rf(x);return RF(add(mul(self.n,x.d),mul(x.n,self.d)),mul(self.d,x.d))
    __radd__=__add__
    def __neg__(self):return RF(tuple(-x for x in self.n),self.d)
    def __sub__(self,x):return self+-as_rf(x)
    def __mul__(self,x):
        x=as_rf(x);return RF(mul(self.n,x.n),mul(self.d,x.d))
    __rmul__=__mul__
    def __truediv__(self,x):
        x=as_rf(x);return RF(mul(self.n,x.d),mul(self.d,x.n))
    def __rtruediv__(self,x):return as_rf(x)/self
    def __eq__(self,x):
        x=as_rf(x);return mul(self.n,x.d)==mul(x.n,self.d)
    def evaluate(self,a):
        f=lambda p:sum(F(v)*a**i for i,v in enumerate(p))
        return f(self.n)/f(self.d)

def as_rf(x):return x if isinstance(x,RF) else RF(x)

def main():
    a=RF((0,1));rates=[RF(1)]*4+[a]
    def weight(s):
        v=RF(1)
        for word in s:
            prefix=RF(0)
            for job in word:
                prefix=prefix+rates[job];v=v/prefix
        return v
    states=[(p[:k],p[k:]) for p in permutations(range(5)) for k in range(6)
            if balanced((p[:k],p[k:]))]
    target_states=[((2,),(0,1,4,3)),((2,),(0,4,1,3))]
    order_A=(0,0,0,1,1);order_B=(0,0,1,1,1)
    D=(a+1)*(a+2)*(a+3)*(a+4)
    expected_As=[1/(6*(a+3)*(a+4)),1/(2*(a+2)*(a+3)*(a+4))]
    records=[]
    for target,expected_A in zip(target_states,expected_As):
        coefficients={o:RF(0) for o in set(orientation(s) for s in states)}
        coefficients[orientation(target)]=-(a+4)*weight(target)
        predecessors=[]
        for s in states:
            for pos,t in enumerate(transitions(s)):
                if t==target:
                    job=(s[0]+s[1])[pos]
                    coefficients[orientation(s)]+=rates[job]*weight(s)
                    predecessors.append({'state':s,'job_served':job,'orientation':orientation(s),
                         'rate_at_a2':str(rates[job].evaluate(F(2))),
                         'weight_at_a2':str(weight(s).evaluate(F(2)))})
        for o,value in coefficients.items():
            expectation=expected_A if o==order_A else (-1/D if o==order_B else RF(0))
            check(value==expectation,'symbolic balance identity mismatch')
        records.append({'state':target,'orientation_A':order_A,'orientation_B':order_B,
                        'predecessors':predecessors})
    ratio1=(a+1)*(a+2)/6;ratio2=(a+1)/2
    check(ratio1-ratio2==(a+1)*(a-1)/6,'ratio incompatibility factorisation failed')
    result={'status':'PASS','parameter':'a>0','class_rates':['1','1','1','1','a'],
      'balance_equations':['x_A/(6(a+3)(a+4)) - x_B/D = 0',
                           'x_A/(2(a+2)(a+3)(a+4)) - x_B/D = 0'],
      'D':'(a+1)(a+2)(a+3)(a+4)',
      'ratio_difference':'(a+1)(a-1)/6',
      'conclusion':'Positive orientation-weighted original product form requires a=1.',
      'method':'Exact polynomial identities by cross multiplication, not parameter sampling.',
      'rows':records}
    Path(__file__).with_name('symbolic_rate_certificate.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:v for k,v in result.items() if k!='rows'},indent=2))

if __name__=='__main__':main()
