"""Compare independently written fixed-position and mutable-list scans.

The mutable scan below is the discovery algorithm copied without its optional
NetworkX atlas dependency. The fixed-position scan is imported from the primary
agent's independent implementation. Only the Python standard library is needed.
"""
from itertools import permutations
import json
from pathlib import Path
from independent_check import transfer

def mutable_move(q, cap, adj, pos=0):
    q = list(q)
    carrier = q.pop(pos)
    scan = pos
    for _ in range(cap):
        while scan < len(q) and not (adj[carrier] >> q[scan]) & 1:
            scan += 1
        if scan == len(q):
            break
        carrier, q[scan] = q[scan], carrier
        scan += 1
    return tuple(q), carrier

def main():
    count=0
    n=5
    adj=[(1<<((i+1)%n))|(1<<((i-1)%n)) for i in range(n)]
    for p in permutations(range(n)):
        for k in range(n+1):
            state=(p[:k],p[k:])
            for source,target in [state,state[::-1]]:
                for pos in range(len(source)):
                    new_q,job=mutable_move(source,2,adj,pos)
                    expected=(new_q,target+(job,))
                    actual=transfer(source,target,pos,cap=2,n=5)
                    if expected != actual:
                        raise RuntimeError(f'Transition mismatch: {state}, position {pos}, {expected}, {actual}')
                    count+=1
    report={'status':'PASS','state_count':720,'event_comparisons':count,
            'swapping_graph':'C5','cap':2,'dependencies':'Python standard library'}
    if count!=3600:
        raise RuntimeError(f'Expected 3600 event comparisons; observed {count}.')
    Path(__file__).with_name('transition_comparison.json').write_text(json.dumps(report,indent=2))
    print(json.dumps(report))

if __name__=='__main__':main()
