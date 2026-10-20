"""The balanced class B_w on C_{2w+1}: closure, irreducibility, law.

Certificate: balanced_class_check.json.  Standard library only; exact;
deterministic.  Run:  python3 -S -O balanced_class_check.py

For w = 2..6 (|B_w| = 180, 1120, 6300, 33264, 168168):
  * every transition from every state of B_w stays in B_w (closure);
  * B_w is strongly connected (forward and backward search from one state);
  * every state has placement height w+1;
  * binding transitions occur only at the head of a full queue, and their
    number is 2Mw/(2w-1), M = (2w+1)!/(w!)^2;
  * |B_w| = (2w+2)!/(w!)^2.
For w = 2, 3, 4 and (mu0, nu0) = (2, 3) the law pi = Bin(n, p) x uniform word,
p = nu0/(mu0+nu0), satisfies global balance exactly at every state.
"""
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations
from math import comb, factorial
from pathlib import Path

HERE = Path(__file__).resolve().parent

def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def dump(payload):
    """Deterministic JSON: one top-level key per line, compact values."""
    return '{\n' + ',\n'.join(json.dumps(k) + ': ' + json.dumps(payload[k], sort_keys=True, separators=(',', ':'))
                              for k in sorted(payload)) + '\n}\n'


def capped_transfer(source, position, adj, w):
    """Serve source[position]; at most w carrier replacements; return
    (new source tuple, departing job).  Direct implementation of the rule."""
    out = list(source)
    carried = source[position]
    loc = position
    for _ in range(w):
        nxt = None
        for j in range(loc + 1, len(source)):
            if adj[carried] >> source[j] & 1:
                nxt = j
                break
        if nxt is None:
            break
        out[nxt] = carried
        carried = source[nxt]
        loc = nxt
    del out[position]
    return tuple(out), carried


def word_height(word, adj):
    depth = {}
    best = 0
    for v in word:
        d = 0
        for u, du in depth.items():
            if adj[v] >> u & 1 and du + 1 > d:
                d = du + 1
        depth[v] = d
        best = max(best, d)
    return best




def balanced_words(w):
    n = 2 * w + 1
    words = []
    for s in range(n):
        for sign in (1, -1):
            long_ = [(s + sign * j) % n for j in range(1, w + 1)]
            short = [(s - sign * j) % n for j in range(1, w)]
            t = (s + sign * (w + 1)) % n
            for pos in combinations(range(2 * w - 1), w):
                ps = set(pos)
                inner = []
                ix = iy = 0
                for j in range(2 * w - 1):
                    if j in ps:
                        inner.append(long_[ix])
                        ix += 1
                    else:
                        inner.append(short[iy])
                        iy += 1
                words.append(tuple([s] + inner + [t]))
    return sorted(set(words))


def uncapped_binds(source, position, adj, w):
    """True if the uncapped scan would make more than w replacements."""
    carried = source[position]
    loc = position
    m = 0
    while True:
        nxt = None
        for j in range(loc + 1, len(source)):
            if adj[carried] >> source[j] & 1:
                nxt = j
                break
        if nxt is None:
            return m > w
        m += 1
        carried = source[nxt]
        loc = nxt


def main():
    rows = []
    for w in range(2, 7):
        n = 2 * w + 1
        adj = [(1 << ((v + 1) % n)) | (1 << ((v - 1) % n)) for v in range(n)]
        words = balanced_words(w)
        M = len(words)
        check(M == factorial(n) // factorial(w) ** 2, 'balanced word count')
        states = [(W[:k], tuple(reversed(W[k:]))) for W in words for k in range(n + 1)]
        check(len(states) == factorial(2 * w + 2) // factorial(w) ** 2, 'class size formula')
        sset = set(states)
        succ = {}
        binding = 0
        binding_at_full_head = 0
        for s in states:
            c, d = s
            check(word_height(c + tuple(reversed(d)), adj) == w + 1, 'height should be w+1')
            out = []
            for p in range(len(c)):
                nc, dep = capped_transfer(c, p, adj, w)
                out.append((nc, d + (dep,)))
                if uncapped_binds(c, p, adj, w):
                    binding += 1
                    binding_at_full_head += (p == 0 and not d)
            for p in range(len(d)):
                nd, dep = capped_transfer(d, p, adj, w)
                out.append((c + (dep,), nd))
                if uncapped_binds(d, p, adj, w):
                    binding += 1
                    binding_at_full_head += (p == 0 and not c)
            for t in out:
                check(t in sset, 'B_w is not closed')
            succ[s] = out
        check(binding == binding_at_full_head, 'binding transition away from a full-queue head')
        check(binding * (2 * w - 1) == 2 * M * w, 'binding count formula')
        pred = {s: [] for s in states}
        for s, ts in succ.items():
            for t in ts:
                pred[t].append(s)
        for graph in (succ, pred):
            seen = {states[0]}
            todo = [states[0]]
            while todo:
                v = todo.pop()
                for u in graph[v]:
                    if u not in seen:
                        seen.add(u)
                        todo.append(u)
            check(len(seen) == len(states), 'B_w is not strongly connected')
        row = {'w': w, 'n': n, 'balanced_words_M': M, 'states': len(states),
               'binding_transitions': binding, 'closed': True, 'strongly_connected': True,
               'height': w + 1}
        if w <= 4:
            mu0, nu0 = F(2), F(3)
            p = nu0 / (mu0 + nu0)
            pi = {s: F(comb(n, len(s[0]))) * p ** len(s[0]) * (1 - p) ** (n - len(s[0])) / M for s in states}
            check(sum(pi.values()) == 1, 'normalisation')
            inflow = {s: F(0) for s in states}
            for s in states:
                c, d = s
                for i, t in enumerate(succ[s]):
                    inflow[t] += pi[s] * (mu0 if i < len(c) else nu0)
            for s in states:
                c, d = s
                check(inflow[s] == pi[s] * (mu0 * len(c) + nu0 * len(d)), 'global balance fails')
            row['homogeneous_law_global_balance_mu0_2_nu0_3'] = True
        rows.append(row)
    payload = {
        'status': 'PASS',
        'rows': rows,
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    (HERE / 'balanced_class_check.json').write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'sizes': [r['states'] for r in rows]}))


if __name__ == '__main__':
    main()
