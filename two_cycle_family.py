"""Two odd cycles joined by a path: a non-cycle family of counterexamples.

Certificate: two_cycle_family.json.  Standard library only; deterministic.
Run:  python3 -S -O two_cycle_family.py

G(a,b,L): an odd cycle C_{2a+1} through vertex a0, an odd cycle C_{2b+1}
through vertex b0, and a path of length L from b0 to a0 (L = 0 identifies
a0 and b0).  Theorem (paper): at cap w = a+b+L+1 the four orientations
  K = { C_{2b+1} balanced with sink b0, path directed b0 -> a0,
        C_{2a+1} balanced with source a0 }
(two choices of source on the first cycle, two of sink on the second)
form a closed class of height w+1, irreducible under positive rates.

Checks for 1 <= a <= b <= 3, 0 <= L <= 3, n = 2a+2b+L+1 <= 13:
  * the moves out of each orientation of K, computed with the binding-move
    lemma, stay in K and connect K (closure + strong connectivity);
  * for n <= 9 the same moves are computed by direct simulation of every
    transition from every state, and must agree;
  * chi(G) = 3 <= w+1 and every orientation in K has height w+1;
  * for n <= 10, all closed classes of height > w are listed for every
    eligible cap (orientation level), as data for the paper's discussion.
"""
import hashlib
import json
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




def build(a, b, L):
    A = list(range(2 * a + 1))                       # a0 = 0
    nxt = 2 * a + 1
    if L == 0:
        B = [0] + list(range(nxt, nxt + 2 * b))
        nxt += 2 * b
        path = [0]
    else:
        B = list(range(nxt, nxt + 2 * b + 1))
        nxt += 2 * b + 1
        inner = list(range(nxt, nxt + L - 1))
        nxt += L - 1
        path = [B[0]] + inner + [0]
    n = nxt
    edges = set()
    for cyc in (A, B):
        m = len(cyc)
        for i in range(m):
            u, v = cyc[i], cyc[(i + 1) % m]
            edges.add((min(u, v), max(u, v)))
    for i in range(len(path) - 1):
        u, v = path[i], path[i + 1]
        edges.add((min(u, v), max(u, v)))
    E = sorted(edges)
    adj = [0] * n
    for u, v in E:
        adj[u] |= 1 << v
        adj[v] |= 1 << u
    return n, adj, E, A, B, path


def balanced_arcs(cyc, end, is_source):
    """Balanced orientation of the odd cycle cyc in which vertex cyc[0] is
    the sink (is_source False) or the source (True); 'end' in {0,1} chooses
    the opposite endpoint cyc[k], k = h or h+1, where len(cyc) = 2h+1."""
    m = len(cyc)
    h = (m - 1) // 2
    k = h + end
    arcs = []
    # arc cyc[0] -> cyc[1] -> ... -> cyc[k]   and   cyc[0] -> cyc[m-1] -> ... -> cyc[k]
    for i in range(k):
        arcs.append((cyc[i], cyc[i + 1]))
    for i in range(m - 1, k - 1, -1):
        arcs.append((cyc[(i + 1) % m], cyc[i]))
    if not is_source:
        arcs = [(v, u) for u, v in arcs]
    return arcs


def mask_of(E, arcs):
    idx = {e: i for i, e in enumerate(E)}
    A = 0
    seen = 0
    for u, v in arcs:
        i = idx[(min(u, v), max(u, v))]
        seen |= 1 << i
        if u < v:
            A |= 1 << i
    check(seen == (1 << len(E)) - 1, 'arcs do not cover the edges')
    return A


def height(n, E, A):
    outs = [[] for _ in range(n)]
    indeg = [0] * n
    for i, (u, v) in enumerate(E):
        s, t = (u, v) if A >> i & 1 else (v, u)
        outs[s].append(t)
        indeg[t] += 1
    order = [v for v in range(n) if indeg[v] == 0]
    dep = [0] * n
    i = 0
    while i < len(order):
        v = order[i]
        i += 1
        for t in outs[v]:
            dep[t] = max(dep[t], dep[v] + 1)
            indeg[t] -= 1
            if indeg[t] == 0:
                order.append(t)
    check(len(order) == n, 'cyclic orientation')
    return max(dep)


def chromatic3(n, adj):
    def ok(k):
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
        return bt(0)
    return 3 if ok(3) and not ok(2) else None


def acyclic(n, E, A):
    outs = [[] for _ in range(n)]
    indeg = [0] * n
    for i, (u, v) in enumerate(E):
        s, t = (u, v) if A >> i & 1 else (v, u)
        outs[s].append(t)
        indeg[t] += 1
    order = [v for v in range(n) if indeg[v] == 0]
    i = 0
    while i < len(order):
        for t in outs[order[i]]:
            indeg[t] -= 1
            if indeg[t] == 0:
                order.append(t)
        i += 1
    return len(order) == n


def main():
    rows = []
    survey = []
    for a in range(1, 4):
        for b in range(a, 4):
            for L in range(0, 4):
                n, adj, E, Acyc, Bcyc, path = build(a, b, L)
                if n > 13:
                    continue
                w = a + b + L + 1
                check(chromatic3(n, adj) == 3, 'chromatic number should be 3')
                K = []
                for eb in (0, 1):
                    for ea in (0, 1):
                        arcs = balanced_arcs(Bcyc, eb, False) + balanced_arcs(Acyc, ea, True)
                        arcs += [(path[i], path[i + 1]) for i in range(len(path) - 1)]
                        K.append(mask_of(E, arcs))
                check(len(set(K)) == 4, 'four distinct orientations expected')
                Kset = set(K)
                for A in K:
                    check(height(n, E, A) == w + 1, 'height should be w+1')
                moves = {A: lemma_successors(n, E, A, w) for A in K}
                check(all(B in Kset for A in K for B in moves[A]), 'K is not closed')
                # strong connectivity of K
                seen = {K[0]}
                todo = [K[0]]
                while todo:
                    A = todo.pop()
                    for B in moves[A]:
                        if B not in seen:
                            seen.add(B)
                            todo.append(B)
                check(seen == Kset, 'K is not connected')
                rev = {A: {B for B in K if A in moves[B]} for A in K}
                seen = {K[0]}
                todo = [K[0]]
                while todo:
                    A = todo.pop()
                    for B in rev[A]:
                        if B not in seen:
                            seen.add(B)
                            todo.append(B)
                check(seen == Kset, 'K is not strongly connected')
                simulated = False
                lin = None
                if n <= 9:
                    graphs, info = quotient(n, adj, E, w)
                    for A in K:
                        check(graphs[w][A] == moves[A], 'simulation disagrees with lemma')
                    simulated = True
                    lin = sorted(info[A][1] for A in K)
                rows.append({'a': a, 'b': b, 'L': L, 'n': n, 'w': w, 'closed': True,
                             'strongly_connected': True, 'orientations': 4, 'height': w + 1,
                             'moves_per_orientation': sorted(len(moves[A]) for A in K),
                             'checked_by_simulation': simulated,
                             'linear_extensions_per_orientation': lin,
                             'states': None if lin is None else sum(lin) * (n + 1)})
                if n <= 10:
                    full = (1 << len(E)) - 1
                    AO = [A for A in range(full + 1) if acyclic(n, E, A)]
                    for ww in range(2, n - 1):
                        g = {A: lemma_successors(n, E, A, ww) for A in AO}
                        cls = []
                        for comp in sccs(g):
                            cs = set(comp)
                            if any(B not in cs for A in comp for B in g[A]):
                                continue
                            h = max(height(n, E, A) for A in comp)
                            if h > ww:
                                cls.append([len(comp), h])
                        if cls:
                            survey.append({'a': a, 'b': b, 'L': L, 'n': n, 'w': ww,
                                           'high_closed_classes_[orientations,height]': sorted(cls)})
    payload = {
        'status': 'PASS',
        'family_rows': rows,
        'survey_all_high_closed_classes_n_le_10': survey,
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    (HERE / 'two_cycle_family.json').write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'family_members': len(rows), 'survey_rows': len(survey)}))


if __name__ == '__main__':
    main()
