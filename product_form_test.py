"""Exact test of placement-order mixtures of original product factors."""
from independent_check import transitions, balanced, orientation, check
from itertools import permutations
from fractions import Fraction as F
from pathlib import Path
import json

def rref(a):
    a=[row[:] for row in a]
    ids=list(range(len(a)));pivots=[];row=0
    for col in range(len(a[0])):
        candidates=[i for i in range(row,len(a)) if a[i][col]]
        if not candidates:continue
        k=min(candidates,key=lambda i:sum(bool(v) for v in a[i]))
        a[row],a[k]=a[k],a[row];ids[row],ids[k]=ids[k],ids[row]
        pivot=a[row][col];a[row]=[v/pivot for v in a[row]]
        for i in range(row+1,len(a)):
            if a[i][col]:
                mul=a[i][col];a[i]=[v-mul*u for v,u in zip(a[i],a[row])]
        pivots.append({'column':col,'original_row':ids[row]});row+=1
        if row==len(a):break
    return row,pivots

def main():
    states=[(p[:k],p[k:]) for p in permutations(range(5)) for k in range(6)
            if balanced((p[:k],p[k:]))]
    si={s:i for i,s in enumerate(states)}
    orders=sorted({orientation(s) for s in states});oi={o:i for i,o in enumerate(orders)}
    all_results=[]
    for rates in [[1]*5,[1,1,1,1,2]]:
        weights=[]
        for s in states:
            value=F(1)
            for word in s:
                prefix=0
                for job in word:
                    prefix+=rates[job];value/=prefix
            weights.append(value)
        mat=[[F(0) for o in orders] for s in states]
        for i,s in enumerate(states):
            col=oi[orientation(s)]
            mat[i][col]-=sum(rates)*weights[i]
            for rate,t in zip([rates[v] for v in s[0]+s[1]],transitions(s)):
                mat[si[t]][col]+=rate*weights[i]
        rank,pivots=rref(mat)
        expected_rank=9 if rates==[1]*5 else 10
        check(rank==expected_rank, 'unexpected exact product-form rank')
        rows=[p['original_row'] for p in pivots]
        all_results.append({'rates':rates,'equation_count':len(mat),'variables':len(orders),
                            'rank':rank,'nullity':len(orders)-rank,
                            'orders':orders,'rank_witness_rows':[
          {'state':states[i],'coefficients':[str(x) for x in mat[i]]} for i in rows]})
    Path(__file__).with_name('product_form_rank.json').write_text(json.dumps(all_results,indent=2)+'\n')
    print(json.dumps([{k:v for k,v in r.items() if k not in ['orders','rank_witness_rows']} for r in all_results],indent=2))

if __name__=='__main__':main()
