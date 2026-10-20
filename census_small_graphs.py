"""Census of capped pass-and-swap closed tandems on all small connected graphs.

Certificate: census_small_graphs.json.  Standard library only; deterministic.
Run:  python3 -S -O census_small_graphs.py

Model: one job per vertex of a simple graph G, cap w at both queues, strictly
positive service rate at every position.  A state is a global word W (a
permutation of the vertices) with a split k (upper queue = W[:k]).

What is computed, for every connected graph on 3..7 vertices (generated here
and counted against OEIS A000088/A001349) and every cap w with
chi(G) <= w+1 <= n-1:
  * every binding transition from every state, by direct simulation of the
    carrier scan, recorded as a move between placement orientations;
  * the closed strongly connected components of this orientation digraph.
Under positive rates all states with one orientation communicate (Lemma
"fixed orientation" of the paper), so a closed component containing an
orientation of height > w is exactly a recurrent class of height > w,
i.e. a violation of Dorsman--Gardner Conjecture 1.

Two cross-checks are also run:
  (a) for all connected graphs on <= 6 vertices and every cap, the closed
      classes of the FULL state graph (no quotient) are computed and their
      high-height members compared with the quotient result;
  (b) for all graphs on <= 6 vertices and every cap, the orientation moves
      predicted by the binding-move lemma of the paper are compared with the
      simulated ones, orientation by orientation.
Failures raise RuntimeError (checks stay active under -O).
"""
import hashlib
import json
import sys
from itertools import permutations
from pathlib import Path

HERE = Path(__file__).resolve().parent
CERT = 'census_small_graphs.json'


def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def dump(payload):
    """Deterministic JSON: one top-level key per line, compact values."""
    return '{\n' + ',\n'.join(json.dumps(k) + ': ' + json.dumps(payload[k], sort_keys=True, separators=(',', ':'))
                              for k in sorted(payload)) + '\n}\n'


# ---------------------------------------------------------------- graphs
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


def all_graphs(nmax):
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
                seen.add(canon(n, adj))
        levels[n] = sorted(seen)
    return levels


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


def chromatic(n, adj):
    for k in range(1, n + 1):
        col = [-1] * n

        def bt(i):
            if i == n:
                return True
            for c in range(k):
                if all(col[u] != c for u in range(i) if adj[i] >> u & 1):
                    col[i] = c
                    if bt(i + 1):
                        return True
            col[i] = -1
            return False
        if bt(0):
            return k
    return n


def edge_list(n, adj):
    return [(u, v) for u in range(n) for v in range(u + 1, n) if adj[u] >> v & 1]


# ---------------------------------------------------------------- dynamics
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


def orientation_bits(word, E):
    pos = {v: i for i, v in enumerate(word)}
    m = 0
    for i, (u, v) in enumerate(E):
        if pos[u] < pos[v]:
            m |= 1 << i
    return m


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


# -------------------------------------------------------------- cross-checks
def full_state_check(n, adj, E, w, info):
    """Closed SCCs of the full state graph; returns sorted list of
    (state count, height) for closed classes containing a high state, and
    asserts that every capped transition's orientation matches the
    quotient's prediction."""
    states = [(p[:k], p[k:]) for p in permutations(range(n)) for k in range(n + 1)]
    idx = {s: i for i, s in enumerate(states)}
    g = {}
    for s in states:
        a, b = s
        nxt = set()
        for p in range(len(a)):
            na, dep = capped_transfer(a, p, adj, w)
            nxt.add(idx[(na, b + (dep,))])
        for p in range(len(b)):
            nb, dep = capped_transfer(b, p, adj, w)
            nxt.add(idx[(a + (dep,), nb)])
        g[idx[s]] = nxt
    word = [s[0] + s[1][::-1] for s in states]
    orient = [orientation_bits(x, E) for x in word]
    res = []
    for comp in sccs(g):
        cs = set(comp)
        if any(t not in cs for s in comp for t in g[s]):
            continue
        hs = [word_height(word[s], adj) for s in comp]
        if max(hs) > w:
            res.append([len(comp), max(hs)])
        # every state of a closed class with one orientation set: check that
        # the class is a union of whole orientation blocks
        oset = {orient[s] for s in comp}
        check(len(comp) == sum(info[A][1] for A in oset) * (n + 1),
              'closed class is not a union of orientation blocks')
    return sorted(res), g, orient


def lemma_successors(n, E, A, w):
    """Binding-move lemma: from orientation A the binding moves are
    flip(x,S) with depth(x)>=w, S a nonempty subset of N+(x) with
    S = N+(x) & down(S+{x}); and the mirror image for the lower queue."""
    outs = [0] * n
    ins = [0] * n
    for i, (u, v) in enumerate(E):
        if A >> i & 1:
            outs[u] |= 1 << v
            ins[v] |= 1 << u
        else:
            outs[v] |= 1 << u
            ins[u] |= 1 << v
    indeg = [bin(ins[v]).count('1') for v in range(n)]
    order = [v for v in range(n) if indeg[v] == 0]
    i = 0
    while i < len(order):
        v = order[i]
        i += 1
        for u in range(n):
            if outs[v] >> u & 1:
                indeg[u] -= 1
                if indeg[u] == 0:
                    order.append(u)
    check(len(order) == n, 'orientation is not acyclic')
    dep = [0] * n
    anc = [0] * n
    for v in order:
        for u in range(n):
            if ins[v] >> u & 1:
                dep[v] = max(dep[v], dep[u] + 1)
                anc[v] |= anc[u] | 1 << u
    up = [0] * n
    desc = [0] * n
    for v in reversed(order):
        for u in range(n):
            if outs[v] >> u & 1:
                up[v] = max(up[v], up[u] + 1)
                desc[v] |= desc[u] | 1 << u
    ebit = {}
    for j, (u, v) in enumerate(E):
        ebit[(u, v)] = ebit[(v, u)] = 1 << j
    res = set()
    for x in range(n):
        for nbset, closure, ok in ((outs[x], anc, dep[x] >= w), (ins[x], desc, up[x] >= w)):
            if not ok or not nbset:
                continue
            nb = [u for u in range(n) if nbset >> u & 1]
            for mask in range(1, 1 << len(nb)):
                S = 0
                for t, u in enumerate(nb):
                    if mask >> t & 1:
                        S |= 1 << u
                D = (1 << x) | closure[x]
                for u in range(n):
                    if S >> u & 1:
                        D |= (1 << u) | closure[u]
                if nbset & D == S:
                    f = 0
                    for u in range(n):
                        if S >> u & 1:
                            f |= ebit[(x, u)]
                    res.add(A ^ f)
    return res


def main():
    nmax = 7
    levels = all_graphs(nmax)
    counts = [len(levels[n]) for n in range(1, nmax + 1)]
    check(counts == [1, 2, 4, 11, 34, 156, 1044], 'graph counts differ from OEIS A000088')
    conn = {n: [c for c in levels[n] if connected(n, decode(n, c))] for n in levels}
    check([len(conn[n]) for n in range(1, nmax + 1)] == [1, 1, 2, 6, 21, 112, 853],
          'connected graph counts differ from OEIS A001349')
    table = {}
    violators = []
    for n in range(3, nmax + 1):
        for code in conn[n]:
            adj = decode(n, code)
            E = edge_list(n, adj)
            chi = chromatic(n, adj)
            ws = [w for w in range(1, n - 1) if chi <= w + 1]
            if not ws:
                continue
            graphs, info = quotient(n, adj, E, max(ws))
            for w in ws:
                key = f'{n},{w}'
                row = table.setdefault(key, {'eligible': 0, 'violating': 0,
                                             'bipartite_eligible': 0, 'bipartite_violating': 0})
                row['eligible'] += 1
                row['bipartite_eligible'] += chi == 2
                bad = high_closed_classes(graphs[w], info, w, n)
                if bad:
                    row['violating'] += 1
                    row['bipartite_violating'] += chi == 2
                    violators.append({'n': n, 'code': code, 'edges': E, 'chi': chi, 'w': w,
                                      'classes_[orientations,states,minheight,maxheight]': bad})
    total = sum(r['violating'] for r in table.values())
    check(total == 1355, f'unexpected number of violating pairs {total}')
    check(all(r['bipartite_violating'] == 0 for r in table.values()),
          'a bipartite violation was found')
    check(all(int(k.split(',')[1]) >= 2 for k, r in table.items() if r['violating']),
          'a violation at cap 1 was found')
    small = sum(r['violating'] for k, r in table.items() if int(k.split(',')[0]) <= 6)
    check(small == 79, 'unexpected number of violating pairs on <= 6 vertices')
    c5 = [v for v in violators if v['n'] == 5 and v['w'] == 2]
    check(len(c5) == 1 and c5[0]['classes_[orientations,states,minheight,maxheight]'] == [[10, 180, 3, 3]],
          'five-cycle class not reproduced')
    # (a) full state graph cross-check, n <= 6, every cap 1..n-1
    full_checked = 0
    for n in range(3, 7):
        for code in conn[n]:
            adj = decode(n, code)
            E = edge_list(n, adj)
            graphs, info = quotient(n, adj, E, n - 1)
            for w in range(1, n):
                fs, g, orient = full_state_check(n, adj, E, w, info)
                qs = sorted([c[1], c[3]] for c in high_closed_classes(graphs[w], info, w, n))
                check(fs == qs, f'full-state and quotient disagree at n={n} code={code} w={w}')
                # transition-level orientation agreement
                realized = {A: set() for A in info}
                for s, ts in g.items():
                    for t in ts:
                        if orient[t] != orient[s]:
                            realized[orient[s]].add(orient[t])
                check(realized == graphs[w], 'quotient moves differ from simulated moves')
                full_checked += 1
    # (b) binding-move lemma versus simulation, all graphs n <= 6
    lemma_checked = 0
    for n in range(2, 7):
        for code in levels[n]:
            adj = decode(n, code)
            E = edge_list(n, adj)
            if not E:
                continue
            graphs, info = quotient(n, adj, E, n - 1)
            for w in range(1, n):
                for A in info:
                    check(graphs[w][A] == lemma_successors(n, E, A, w),
                          f'binding-move lemma fails at n={n} code={code} w={w}')
                    lemma_checked += 1
    payload = {
        'status': 'PASS',
        'model': 'one job per vertex, cap w at both queues, positive rate at every position',
        'graph_counts_all_n1_to_7': counts,
        'connected_counts_n1_to_7': [len(conn[n]) for n in range(1, nmax + 1)],
        'table_by_n_w': table,
        'violating_pairs_total': total,
        'violating_pairs_n_le_6': small,
        'violators': violators,
        'full_state_crosscheck_graph_cap_pairs': full_checked,
        'binding_move_lemma_orientation_checks': lemma_checked,
        'code_convention': 'canonical upper-triangle bit code: bit index runs over pairs (i,j), i<j, lexicographically',
    }
    src = Path(__file__).read_bytes()
    payload['input_sha256'] = hashlib.sha256(src).hexdigest()
    payload['verifier'] = Path(__file__).name
    (HERE / CERT).write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'violating_pairs_total': total,
                      'full_state_crosscheck_graph_cap_pairs': full_checked,
                      'binding_move_lemma_orientation_checks': lemma_checked}))


if __name__ == '__main__':
    main()
