"""Independent finite-state reconstruction by the primary agent.

No code is imported from the discovery implementation.  Jobs use labels 0,...,4.
"""
from itertools import permutations
from collections import Counter, deque
from fractions import Fraction
import json
from math import factorial
from pathlib import Path

def check(condition, message):
    if not condition:
        raise RuntimeError(message)

def transfer(source, target, position, cap=2, n=5):
    out=list(source)
    carried=source[position]
    loc=position
    for count in range(cap):
        candidates=[j for j in range(loc+1,len(source))
                    if (source[j]-carried)%n in (1,n-1)]
        if not candidates:
            break
        nxt=min(candidates)
        out[nxt]=carried
        carried=source[nxt]
        loc=nxt
    del out[position]
    return tuple(out),target+(carried,)

def transitions(state):
    a,b=state
    ans=[]
    for p in range(len(a)):
        ans.append(transfer(a,b,p))
    for p in range(len(b)):
        bp,ap=transfer(b,a,p)
        ans.append((ap,bp))
    return ans

def orientation(s):
    word=s[0]+s[1][::-1]
    rank={v:i for i,v in enumerate(word)}
    return tuple(int(rank[v]<rank[(v+1)%5]) for v in range(5))

def balanced(s):
    o=orientation(s)
    return sum(o) in (2,3) and sum(o[i]!=o[(i+1)%5] for i in range(5))==2

def placement_height(s):
    word=s[0]+s[1][::-1]
    distance={}
    for vertex in word:
        predecessors=[v for v in ((vertex-1)%5,(vertex+1)%5) if v in distance]
        distance[vertex]=max((distance[v]+1 for v in predecessors), default=0)
    return max(distance.values())

def canonical(s):
    versions=[]
    for sign in (-1,1):
        for shift in range(5):
            a,b=(tuple((sign*v+shift)%5 for v in row) for row in s)
            versions.extend([(a,b),(b,a)])
    return min(versions)

def reachable(start, graph, allowed=None):
    seen={start};todo=[start]
    while todo:
        v=todo.pop()
        for u in graph[v]:
            if u not in seen and (allowed is None or u in allowed):
                seen.add(u);todo.append(u)
    return seen

def components(graph):
    reverse={s:[] for s in graph}
    for s,ts in graph.items():
        for t in ts:reverse[t].append(s)
    remaining=set(graph);result=[]
    while remaining:
        seed=min(remaining)
        component=reachable(seed,graph,remaining)&reachable(seed,reverse,remaining)
        result.append(component);remaining-=component
    return result

def solve_fraction(matrix,rhs):
    n=len(rhs)
    a=[[Fraction(v) for v in row]+[Fraction(rhs[i])] for i,row in enumerate(matrix)]
    for col in range(n):
        pivot=next(i for i in range(col,n) if a[i][col])
        a[col],a[pivot]=a[pivot],a[col]
        scale=a[col][col]
        a[col]=[v/scale for v in a[col]]
        for i in range(n):
            if i!=col:
                scale=a[i][col]
                a[i]=[v-scale*w for v,w in zip(a[i],a[col])]
    return [a[i][-1] for i in range(n)]

def main():
    states=[(p[:k],p[k:]) for p in permutations(range(5)) for k in range(6)]
    chosen=[s for s in states if balanced(s)]
    selected=set(chosen)
    check(len(states)==720, 'incorrect full state count')
    check(len(chosen)==180, 'incorrect balanced class size')
    adj={s:transitions(s) for s in states}
    check(all(t in selected for s in chosen for t in adj[s]), 'balanced class is not closed')
    check(len(components({s:adj[s] for s in chosen}))==1, 'balanced class is not strongly connected')
    sccs=components(adj)
    closed=[c for c in sccs if all(t in c for s in c for t in adj[s])]
    check(sorted(map(len,closed))==[48]*10+[180], 'unexpected closed SCC sizes')
    check(len(states)-sum(map(len,closed))==60, 'unexpected transient state count')
    heights=Counter(placement_height(s) for s in states)
    check(heights=={2:480,3:180,4:60}, 'unexpected placement-height counts')
    check(all(placement_height(s)==3 for s in chosen), 'balanced class has unexpected height')
    reps=sorted({canonical(s) for s in chosen})
    ri={s:i for i,s in enumerate(reps)}
    sizes=Counter(canonical(s) for s in chosen)
    check(len(reps)==9 and set(sizes.values())=={20}, 'unexpected symmetry orbits')
    counts=[[0]*len(reps) for _ in reps]
    for i,s in enumerate(reps):
        for t in adj[s]: counts[i][ri[canonical(t)]]+=1
    for s in chosen:
        row=Counter(ri[canonical(t)] for t in adj[s])
        check([row[i] for i in range(9)]==counts[ri[canonical(s)]], 'quotient is not lumpable')
    q=[[counts[i][j]-5*(i==j) for j in range(9)] for i in range(9)]
    equations=[[q[i][j] for i in range(9)] for j in range(8)]+[[1]*9]
    stationary=solve_fraction(equations,[0]*8+[1])
    check(all(sum(stationary[i]*q[i][j] for i in range(9))==0 for j in range(9)), 'stationary residual is nonzero')
    check(sum(stationary)==1 and all(v>0 for v in stationary), 'stationary probabilities do not normalise positively')
    # Published original product form for unit-rate service of every job.
    original={s:Fraction(1,factorial(len(s[0]))*factorial(len(s[1]))) for s in chosen}
    residual={s:sum(original[t]*adj[t].count(s) for t in chosen)-5*original[s] for s in chosen}
    witness=next((s for s in chosen if residual[s]),None)
    check(witness is None, 'homogeneous product weights fail global balance')
    normaliser=sum(original.values())
    check(all(original[s]/normaliser==stationary[ri[canonical(s)]]/sizes[canonical(s)] for s in chosen),
          'closed-form and independently solved stationary probabilities disagree')
    rates=[1,1,1,1,2]
    def heterogeneous_weight(s):
        result=Fraction(1)
        for word in s:
            prefix=0
            for v in word:
                prefix+=rates[v];result/=prefix
        return result
    heterogeneous={s:heterogeneous_weight(s) for s in chosen}
    incoming=Counter()
    for s in chosen:
        for rate,t in zip([rates[v] for v in s[0]+s[1]],adj[s]):
            incoming[t]+=rate*heterogeneous[s]
    hetero_residual={s:incoming[s]-sum(rates)*heterogeneous[s] for s in chosen}
    hwitness=next((s for s in chosen if hetero_residual[s]),None)
    check(hwitness is not None, 'expected heterogeneous obstruction not found')
    result={'status':'PASS','arithmetic':'Python integers and fractions; acceptance checks remain active under -O',
      'state_count':len(states),'closed_scc_sizes':sorted(map(len,closed)),
      'transient_state_count':len(states)-sum(map(len,closed)),
      'placement_height_counts':dict(sorted(heights.items())),
      'bad_class_size':len(chosen),'orbits':[{'index':i,'representative':s,'orbit_size':sizes[s],
      'orbit_probability':str(stationary[i]),'per_state_probability':str(stationary[i]/sizes[s]),
      'transition_counts':counts[i]} for i,s in enumerate(reps)],
      'equal_rate_product_form_valid':witness is None,
      'heterogeneous_class_rates':rates,
      'heterogeneous_product_form_residual_witness':None if hwitness is None else
      {'state':hwitness,'unnormalised_weight':str(heterogeneous[hwitness]),
       'balance_residual':str(hetero_residual[hwitness])}}
    Path(__file__).with_name('independent_c5.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__': main()
