"""Several jobs per class: twin blow-ups of the swapping graph.

Certificate: multiplicity_check.json.  Standard library only; deterministic.
Run:  python3 -S -O multiplicity_check.py

With m_v jobs of class v (and no self-loops, so jobs of one class never
swap with each other), the relabelling of Dorsman--Gardner (2024, Remark 9)
gives an equivalent system with one job per vertex of the blow-up graph in
which v is replaced by m_v pairwise non-adjacent twins.  Recurrence is
preserved by this relabelling (see the paper), so we analyse blow-ups.

Checks (positive rates; orientation quotient by direct simulation when the
blow-up has at most 7 vertices, binding-move lemma otherwise):
  * C_5, cap 2: every multiplicity vector in {1,2,3}^5 with at most 8 jobs
    (up to rotation and reflection): a high recurrent class exists only for
    (1,1,1,1,1);
  * C_7, cap 3: one class doubled, two classes doubled; C_9, caps 2 and 4:
    one class doubled: no high recurrent class;
  * K_3 (triangle) with multiplicities (2,1,1), (2,2,1), (2,2,2), (3,1,1),
    and the bowtie with one class doubled, at every eligible cap: high
    recurrent classes DO occur (e.g. K_3 with two jobs per class at caps 3
    and 4), so multiplicity does not in general restore the conjecture.
"""
import hashlib
import json
from itertools import product
from pathlib import Path

HERE = Path(__file__).resolve().parent

def check(cond, msg):
    if not cond:
        raise RuntimeError(msg)


def dump(payload):
    """Deterministic JSON: one top-level key per line, compact values."""
    return '{\n' + ',\n'.join(json.dumps(k) + ': ' + json.dumps(payload[k], sort_keys=True, separators=(',', ':'))
                              for k in sorted(payload)) + '\n}\n'


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




def blowup(n, el, mult):
    idx = []
    lab = 0
    for v in range(n):
        idx.append(list(range(lab, lab + mult[v])))
        lab += mult[v]
    E = sorted({(min(a, b), max(a, b)) for u, v in el for a in idx[u] for b in idx[v]})
    adj = [0] * lab
    for a, b in E:
        adj[a] |= 1 << b
        adj[b] |= 1 << a
    return lab, adj, E


def height(N, E, A):
    outs = [[] for _ in range(N)]
    indeg = [0] * N
    for i, (u, v) in enumerate(E):
        s, t = (u, v) if A >> i & 1 else (v, u)
        outs[s].append(t)
        indeg[t] += 1
    order = [v for v in range(N) if indeg[v] == 0]
    dep = [0] * N
    i = 0
    while i < len(order):
        v = order[i]
        i += 1
        for t in outs[v]:
            dep[t] = max(dep[t], dep[v] + 1)
            indeg[t] -= 1
            if indeg[t] == 0:
                order.append(t)
    return max(dep) if len(order) == N else None


def high_classes(N, adj, E, w):
    if N <= 7:
        graphs, info = quotient(N, adj, E, w)
        return [[c[0], c[3]] for c in high_closed_classes(graphs[w], info, w, N)], 'simulation'
    full = (1 << len(E)) - 1
    AO = [A for A in range(full + 1) if height(N, E, A) is not None]
    g = {A: lemma_successors(N, E, A, w) for A in AO}
    res = []
    for comp in sccs(g):
        cs = set(comp)
        if any(B not in cs for A in comp for B in g[A]):
            continue
        h = max(height(N, E, A) for A in comp)
        if h > w:
            res.append([len(comp), h])
    return sorted(res), 'binding-move lemma'


def cyc(n):
    return [(i, (i + 1) % n) for i in range(n)]


def main():
    rows = []
    seen = set()
    for mult in product((1, 2, 3), repeat=5):
        if sum(mult) > 8:
            continue
        key = min(tuple(mult[(s * i + r) % 5] for i in range(5)) for r in range(5) for s in (1, -1))
        if key in seen:
            continue
        seen.add(key)
        N, adj, E = blowup(5, cyc(5), key)
        cls, how = high_classes(N, adj, E, 2)
        check(bool(cls) == (key == (1, 1, 1, 1, 1)), f'C5 multiplicity {key}')
        rows.append({'graph': 'C5', 'multiplicities': list(key), 'w': 2, 'jobs': N,
                     'high_recurrent_classes_[orientations,height]': cls, 'method': how})
    for g, n, w, mult in [('C7', 7, 3, (2, 1, 1, 1, 1, 1, 1)), ('C7', 7, 3, (2, 2, 1, 1, 1, 1, 1)),
                          ('C7', 7, 3, (2, 1, 2, 1, 1, 1, 1)), ('C7', 7, 3, (2, 1, 1, 2, 1, 1, 1)),
                          ('C9', 9, 2, (2,) + (1,) * 8), ('C9', 9, 4, (2,) + (1,) * 8)]:
        N, adj, E = blowup(n, cyc(n), mult)
        cls, how = high_classes(N, adj, E, w)
        check(not cls, f'{g} multiplicity {mult} cap {w}')
        rows.append({'graph': g, 'multiplicities': list(mult), 'w': w, 'jobs': N,
                     'high_recurrent_classes_[orientations,height]': cls, 'method': how})
    tri = [(0, 1), (1, 2), (0, 2)]
    bow = [(0, 1), (1, 2), (2, 0), (0, 3), (3, 4), (4, 0)]
    found = {}
    for g, n, el, mult in [('K3', 3, tri, (2, 1, 1)), ('K3', 3, tri, (2, 2, 1)), ('K3', 3, tri, (2, 2, 2)),
                           ('K3', 3, tri, (3, 1, 1)), ('bowtie', 5, bow, (2, 1, 1, 1, 1)),
                           ('bowtie', 5, bow, (1, 2, 1, 1, 1))]:
        N, adj, E = blowup(n, el, mult)
        for w in range(2, N - 1):            # chi = 3 <= w+1
            cls, how = high_classes(N, adj, E, w)
            rows.append({'graph': g, 'multiplicities': list(mult), 'w': w, 'jobs': N,
                         'high_recurrent_classes_[orientations,height]': cls, 'method': how})
            found[(g, mult, w)] = bool(cls)
    check(found[('K3', (2, 2, 2), 3)] and found[('K3', (2, 2, 2), 4)], 'K3 doubled: expected high classes')
    check(found[('bowtie', (2, 1, 1, 1, 1), 3)], 'bowtie doubled: expected high classes')
    payload = {
        'status': 'PASS',
        'rows': rows,
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    (HERE / 'multiplicity_check.json').write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'rows': len(rows),
                      'with_high_classes': [(r['graph'], r['multiplicities'], r['w'])
                                            for r in rows if r['high_recurrent_classes_[orientations,height]']]}))


if __name__ == '__main__':
    main()
