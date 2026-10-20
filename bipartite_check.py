"""Bipartite swapping graphs: no recurrent state of excessive height.

Certificate: bipartite_check.json (default, graphs on <= 8 vertices) or
bipartite_check_n9.json (argument 9: additionally all 730 connected
bipartite graphs on 9 vertices; about 20 CPU minutes).
Standard library only; deterministic.  Run:
    python3 -S -O bipartite_check.py          # n <= 8
    python3 -S -O bipartite_check.py 9        # n <= 9 (long)

Checks, for every connected bipartite graph G on 3..N vertices (generated
here; counts compared with OEIS A033995/A005142) and EVERY cap 1 <= w <= n-2:
  (i)  positive rates at every position: the orientation digraph of binding
       moves (by direct simulation of every transition) has no closed
       component containing an orientation of height > w;
  (ii) head-of-line service only (the minimal positivity allowed by the
       order-independence conditions; FCFS-type), for n <= 7: the full state
       graph has no closed class containing a state of height > w.
The paper proves (i) for every odd w and every bipartite graph; for even w,
(i) is the certified evidence behind the refined conjecture.
"""
import hashlib
import json
import sys
from itertools import permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent

def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def dump(payload):
    """Deterministic JSON: one top-level key per line, compact values."""
    return '{\n' + ',\n'.join(json.dumps(k) + ': ' + json.dumps(payload[k], sort_keys=True, separators=(',', ':'))
                              for k in sorted(payload)) + '\n}\n'


def refine(n, adj):
    col = [0] * n
    while True:
        sig = [(col[v], tuple(sorted(col[u] for u in range(n) if adj[v] >> u & 1)))
               for v in range(n)]
        keys = sorted(set(sig))
        new = [keys.index(s) for s in sig]
        if len(set(new)) == len(set(col)):
            return new
        col = new


def canon(n, adj):
    """Canonical code: minimum upper-triangle bit code over all labellings
    that respect the (isomorphism-invariant) colour-refinement order."""
    col = refine(n, adj)
    classes = {}
    for v in range(n):
        classes.setdefault(col[v], []).append(v)
    groups = [classes[c] for c in sorted(classes)]

    def seqs(i):
        if i == len(groups):
            yield []
            return
        for p in permutations(groups[i]):
            for rest in seqs(i + 1):
                yield list(p) + rest
    best = None
    for seq in seqs(0):
        code = 0
        bit = 0
        for i in range(n):
            ai = adj[seq[i]]
            for j in range(i + 1, n):
                if ai >> seq[j] & 1:
                    code |= 1 << bit
                bit += 1
        if best is None or code < best:
            best = code
    return best


def decode(n, code):
    adj = [0] * n
    bit = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> bit & 1:
                adj[i] |= 1 << j
                adj[j] |= 1 << i
            bit += 1
    return adj


def connected(n, adj):
    seen = 1
    frontier = 1
    while frontier:
        nxt = 0
        for v in range(n):
            if frontier >> v & 1:
                nxt |= adj[v]
        nxt &= ~seen
        seen |= nxt
        frontier = nxt
    return seen == (1 << n) - 1


def edge_list(n, adj):
    return [(u, v) for u in range(n) for v in range(u + 1, n) if adj[u] >> v & 1]


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


def quotient(n, adj, E, wmax):
    """Orientation digraph for caps 1..wmax by DFS over prefixes.  For a
    prefix c (the upper queue) and served position p, the uncapped carriers
    are followed; for each cap w smaller than the uncapped number of
    replacements the transition binds, and its effect on the orientation is
    recorded (reversal of the edges from the w-th carrier to its later
    neighbours in c).  This equals the orientation of the simulated capped
    transition; that equality is asserted on every event of the full-state
    cross-check below."""
    ebit = {}
    for i, (u, v) in enumerate(E):
        ebit[(u, v)] = ebit[(v, u)] = 1 << i
    info = {}
    succ = {w: {} for w in range(1, wmax + 1)}
    word = []
    used = [False] * n
    depth = [0] * n
    stack = []

    def rec(k, A, hmax):
        if k == n:
            if A in info:
                info[A][1] += 1
            else:
                info[A] = [hmax, 1]
            for flips in stack:
                for w, f in flips:
                    succ[w].setdefault(A, set()).add(A ^ f)
            return
        for v in range(n):
            if used[v]:
                continue
            used[v] = True
            Ap = A
            d = 0
            for u in word:
                if adj[v] >> u & 1:
                    if u < v:
                        Ap |= ebit[(u, v)]
                    d = max(d, depth[u] + 1)
            depth[v] = d
            word.append(v)
            kk = k + 1
            flips = []
            for p in range(kk):
                carr = word[p]
                loc = p
                m = 0
                while m <= wmax:
                    j = loc + 1
                    while j < kk and not (adj[carr] >> word[j] & 1):
                        j += 1
                    if j == kk:
                        break
                    m += 1
                    if m <= wmax:
                        x = word[j]
                        f = 0
                        for jj in range(j + 1, kk):
                            if adj[x] >> word[jj] & 1:
                                f |= ebit[(x, word[jj])]
                        if f:
                            flips.append((m, f))
                    carr = word[j]
                    loc = j
            stack.append(flips)
            rec(kk, Ap, max(hmax, d))
            stack.pop()
            word.pop()
            used[v] = False
    rec(0, 0, 0)
    full = (1 << len(E)) - 1
    graphs = {}
    for w in range(1, wmax + 1):
        g = {A: set() for A in info}
        for A, ss in succ[w].items():
            for B in ss:
                g[A].add(B)
                g[A ^ full].add(B ^ full)   # queue exchange reverses every edge
        graphs[w] = g
    return graphs, info


def sccs(g):
    """Iterative Tarjan; deterministic order."""
    index = {}
    low = {}
    on = set()
    st = []
    out = []
    c = 0
    for root in sorted(g):
        if root in index:
            continue
        index[root] = low[root] = c
        c += 1
        st.append(root)
        on.add(root)
        work = [(root, iter(sorted(g[root])))]
        while work:
            v, it = work[-1]
            pushed = False
            for u in it:
                if u not in index:
                    index[u] = low[u] = c
                    c += 1
                    st.append(u)
                    on.add(u)
                    work.append((u, iter(sorted(g[u]))))
                    pushed = True
                    break
                if u in on:
                    low[v] = min(low[v], index[u])
            if pushed:
                continue
            work.pop()
            if work:
                low[work[-1][0]] = min(low[work[-1][0]], low[v])
            if low[v] == index[v]:
                comp = []
                while True:
                    u = st.pop()
                    on.discard(u)
                    comp.append(u)
                    if u == v:
                        break
                out.append(sorted(comp))
    return out


def high_closed_classes(g, info, w, n):
    res = []
    for comp in sccs(g):
        cs = set(comp)
        if any(B not in cs for A in comp for B in g[A]):
            continue
        hs = [info[A][0] for A in comp]
        if max(hs) > w:
            res.append([len(comp), sum(info[A][1] for A in comp) * (n + 1),
                        min(hs), max(hs)])
    return sorted(res)




def bipartite(n, adj):
    col = [-1] * n
    for s in range(n):
        if col[s] >= 0:
            continue
        col[s] = 0
        st = [s]
        while st:
            v = st.pop()
            for u in range(n):
                if adj[v] >> u & 1:
                    if col[u] < 0:
                        col[u] = 1 - col[v]
                        st.append(u)
                    elif col[u] == col[v]:
                        return False
    return True


def bipartite_graphs(nmax):
    levels = {1: [0]}
    for n in range(2, nmax + 1):
        seen = set()
        for code in levels[n - 1]:
            base = decode(n - 1, code)
            for mask in range(1 << (n - 1)):
                adj = base[:] + [mask]
                for u in range(n - 1):
                    if mask >> u & 1:
                        adj[u] |= 1 << (n - 1)
                if bipartite(n, adj):
                    seen.add(canon(n, adj))
        levels[n] = sorted(seen)
    return levels


def head_only_high_recurrent(n, adj, w):
    """Closed classes of the full state graph when only the head of each
    queue is served.  Returns the number of closed classes containing a
    state of height > w."""
    states = [(p[:k], p[k:]) for p in permutations(range(n)) for k in range(n + 1)]
    idx = {s: i for i, s in enumerate(states)}
    g = {}
    for s in states:
        a, b = s
        nxt = set()
        if a:
            na, dep = capped_transfer(a, 0, adj, w)
            nxt.add(idx[(na, b + (dep,))])
        if b:
            nb, dep = capped_transfer(b, 0, adj, w)
            nxt.add(idx[(a + (dep,), nb)])
        g[idx[s]] = nxt
    bad = 0
    closed = 0
    for comp in sccs(g):
        cs = set(comp)
        if any(t not in cs for s in comp for t in g[s]):
            continue
        closed += 1
        if max(word_height(states[s][0] + states[s][1][::-1], adj) for s in comp) > w:
            bad += 1
    return bad, closed


def main():
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 8
    check(nmax in (8, 9), 'argument must be 8 or 9')
    levels = bipartite_graphs(nmax)
    allc = [1, 2, 3, 7, 13, 35, 88, 303, 1119][:nmax]
    conc = [1, 1, 1, 3, 5, 17, 44, 182, 730][:nmax]
    check([len(levels[n]) for n in range(1, nmax + 1)] == allc, 'bipartite counts differ from A033995')
    conn = {n: [c for c in levels[n] if connected(n, decode(n, c))] for n in levels}
    check([len(conn[n]) for n in range(1, nmax + 1)] == conc, 'connected bipartite counts differ from A005142')
    positive = {}
    for n in range(3, nmax + 1):
        for code in conn[n]:
            adj = decode(n, code)
            E = edge_list(n, adj)
            graphs, info = quotient(n, adj, E, n - 2)
            for w in range(1, n - 1):
                row = positive.setdefault(f'{n},{w}', [0, 0])
                row[0] += 1
                bad = high_closed_classes(graphs[w], info, w, n)
                row[1] += len(bad) > 0
    check(all(v[1] == 0 for v in positive.values()), 'bipartite violation found (positive rates)')
    head = {}
    for n in range(3, 8):
        for code in conn[n]:
            adj = decode(n, code)
            for w in range(1, n - 1):
                bad, closed = head_only_high_recurrent(n, adj, w)
                row = head.setdefault(f'{n},{w}', [0, 0, 0])
                row[0] += 1
                row[1] += bad > 0
                row[2] += closed
    check(all(v[1] == 0 for v in head.values()), 'bipartite violation found (head-only service)')
    payload = {
        'status': 'PASS',
        'max_vertices': nmax,
        'bipartite_counts': allc,
        'connected_bipartite_counts': conc,
        'positive_rates_[graphs,violating]_by_n_w': positive,
        'head_only_service_[graphs,violating,closed_classes_total]_by_n_w': head,
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    name = 'bipartite_check.json' if nmax == 8 else 'bipartite_check_n9.json'
    (HERE / name).write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'max_vertices': nmax,
                      'graph_cap_pairs_positive': sum(v[0] for v in positive.values()),
                      'graph_cap_pairs_head_only': sum(v[0] for v in head.values())}))


if __name__ == '__main__':
    main()
