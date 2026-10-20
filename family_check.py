"""Exact odd-cycle head-bijection checks; Python standard library only.

This is a finite check of the paper's general proof, not a proof by enumeration.
Run: python family_check.py [maximum_cap], default maximum_cap=7.
"""
from itertools import combinations
from math import comb, factorial
from pathlib import Path
import json
import sys

def balanced_words(w):
    n=2*w+1
    words=set()
    for source in range(n):
        for sign in (-1,1):
            long=[(source+sign*j)%n for j in range(1,w+1)]
            short=[(source-sign*j)%n for j in range(1,w)]
            sink=(source+sign*(w+1))%n
            for positions in combinations(range(2*w-1),w):
                inner=[]
                ix=iy=0
                chosen=set(positions)
                for j in range(2*w-1):
                    if j in chosen:
                        inner.append(long[ix]);ix+=1
                    else:
                        inner.append(short[iy]);iy+=1
                words.add((source,*inner,sink))
    return words

def head_scan(word,w):
    n=len(word)
    q=list(word)
    carried=q.pop(0)
    scan=0
    for _ in range(w):
        while scan<len(q) and (q[scan]-carried)%n not in (1,n-1):
            scan+=1
        if scan==len(q):break
        carried,q[scan]=q[scan],carried
        scan+=1
    return tuple(q)+(carried,)

def main():
    maximum=int(sys.argv[1]) if len(sys.argv)>1 else 7
    if maximum < 2:
        raise RuntimeError('At least cap 2 must be checked.')
    results=[]
    for w in range(2,maximum+1):
        n=2*w+1
        words=balanced_words(w)
        images=[head_scan(word,w) for word in words]
        expected=2*n*comb(2*w-1,w-1)
        if not len(words)==expected==factorial(n)//factorial(w)**2:
            raise RuntimeError(f'Balanced-word count disagrees at cap {w}.')
        if not (len(set(images))==len(images) and set(images)==words):
            raise RuntimeError(f'Capped head map is not a bijection at cap {w}.')
        results.append({'w':w,'n':n,'full_queue_words':len(words),
                        'total_closed_states':len(words)*(n+1),
                        'capped_head_bijection':True})
    path=Path(__file__).with_name('family_head_bijection.json')
    path.write_text(json.dumps(results,indent=2))
    print(json.dumps({'status':'PASS','caps_checked':len(results),
                      'largest_cap':maximum,'output':path.name}))

if __name__=='__main__':main()
