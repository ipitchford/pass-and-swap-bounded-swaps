"""Cycles: which caps admit a recurrent class of excessive height.

Certificate: cycle_classification.json.  Standard library only; deterministic.
Run:  python3 -S -O cycle_classification.py

Claim checked (Theorem "cycles" of the paper, proved there for all n, w):
with one job per vertex of C_n, cap w (1 <= w <= n-2) and positive rates,
a recurrent state of height > w exists iff n = 1 (mod 2w); in that case there
is exactly one such class; its orientations are the 2n orientations whose
runs (maximal directed arcs) have lengths w+1, w, w, ..., w (2r runs,
n = 2rw+1), and its size is 2n(n+1)e, e = number of linear extensions.

  (a) 3 <= n <= 9: by direct simulation of every transition (orientation
      quotient, as in census_small_graphs.py); also the run rule for binding
      moves (flip edge j of a run of length l >= w+1 iff j >= w or
      j <= l-w-1) is compared with the simulated moves for every orientation;
  (b) 10 <= n <= 14: orientation level, using the run rule.
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




def cycle(n):
    adj = [(1 << ((v + 1) % n)) | (1 << ((v - 1) % n)) for v in range(n)]
    E = sorted((min(v, (v + 1) % n), max(v, (v + 1) % n)) for v in range(n))
    return adj, E


def directions(n, E, A):
    """dir[i] = +1 if edge {i, i+1} points i -> i+1 (clockwise), else -1."""
    d = [0] * n
    for b, (u, v) in enumerate(E):
        fwd = bool(A >> b & 1)          # u -> v
        if v == u + 1:
            d[u] = 1 if fwd else -1
        else:                            # (0, n-1): edge n-1 -- 0
            d[n - 1] = -1 if fwd else 1
    return d


def from_directions(n, E, d):
    A = 0
    for b, (u, v) in enumerate(E):
        if v == u + 1:
            fwd = d[u] == 1
        else:
            fwd = d[n - 1] == -1
        if fwd:
            A |= 1 << b
    return A


def runs(d):
    """Maximal runs of equal direction, cyclically: list of (start edge, length, sign)."""
    n = len(d)
    if all(x == d[0] for x in d):
        return None                      # directed cycle (not acyclic)
    s = 0
    while d[s - 1] == d[s]:
        s -= 1
    s %= n
    out = []
    i = 0
    while i < n:
        j = i
        while j + 1 < n and d[(s + j + 1) % n] == d[(s + i) % n]:
            j += 1
        out.append(((s + i) % n, j - i + 1, d[(s + i) % n]))
        i = j + 1
    return out


def run_rule_moves(n, E, A, w):
    d = directions(n, E, A)
    res = set()
    for start, length, sign in runs(d):
        if length < w + 1:
            continue
        for j in range(length):          # j-th edge along the run's direction
            if j >= w or j <= length - w - 1:
                e = (start + j) % n if sign == 1 else (start + length - 1 - j) % n
                d2 = d[:]
                d2[e] = -d2[e]
                res.add(from_directions(n, E, d2))
    return res


def height_from_runs(d):
    return max(length for _, length, _ in runs(d))


def lin_ext(n, arcs):
    pred = [0] * n
    for a, b in arcs:
        pred[b] |= 1 << a
    f = {0: 1}
    for mask in range(1 << n):
        if mask not in f:
            continue
        for v in range(n):
            if not mask >> v & 1 and pred[v] & mask == pred[v]:
                f[mask | 1 << v] = f.get(mask | 1 << v, 0) + f[mask]
    return f[(1 << n) - 1]


def arcs_of(n, d):
    return [(i, (i + 1) % n) if d[i] == 1 else ((i + 1) % n, i) for i in range(n)]


def expected_pattern(n, w):
    if (n - 1) % (2 * w) or n < 2 * w + 1:
        return None
    r = (n - 1) // (2 * w)
    return sorted([w + 1] + [w] * (2 * r - 1))


def main():
    rows = []
    for n in range(3, 15):
        adj, E = cycle(n)
        full = (1 << n) - 1
        if n <= 9:
            graphs, info = quotient(n, adj, E, n - 2)
            orients = sorted(info)
            # run rule == simulated moves
            for w in range(1, n - 1):
                for A in orients:
                    check(graphs[w][A] == run_rule_moves(n, E, A, w),
                          f'run rule fails at n={n} w={w}')
            method = 'simulation of every transition'
        else:
            orients = [A for A in range(full + 1) if runs(directions(n, E, A)) is not None]
            info = None
            method = 'orientation level (run rule)'
        check(len(orients) == 2 ** n - 2, 'acyclic orientation count')
        for w in range(1, n - 1):
            if n <= 9:
                g = graphs[w]
            else:
                g = {A: run_rule_moves(n, E, A, w) for A in orients}
            hts = {A: height_from_runs(directions(n, E, A)) for A in orients}
            bad = []
            for comp in sccs(g):
                cs = set(comp)
                if any(B not in cs for A in comp for B in g[A]):
                    continue
                if max(hts[A] for A in comp) > w:
                    bad.append(comp)
            pat = expected_pattern(n, w)
            if pat is None:
                check(not bad, f'unexpected high class n={n} w={w}')
                rows.append({'n': n, 'w': w, 'high_recurrent_classes': 0, 'method': method})
                continue
            check(len(bad) == 1, f'expected exactly one high class n={n} w={w}')
            comp = bad[0]
            check(len(comp) == 2 * n, 'class should have 2n orientations')
            ext = set()
            for A in comp:
                d = directions(n, E, A)
                check(sorted(L for _, L, _ in runs(d)) == pat, 'run pattern mismatch')
                ext.add(lin_ext(n, arcs_of(n, d)))
            check(len(ext) == 1, 'linear extension counts differ within class')
            e = ext.pop()
            states = 2 * n * (n + 1) * e
            if info is not None:
                check(states == sum(info[A][1] for A in comp) * (n + 1), 'state count mismatch')
            rows.append({'n': n, 'w': w, 'high_recurrent_classes': 1, 'orientations': len(comp),
                         'run_lengths': pat, 'linear_extensions_per_orientation': e,
                         'states': states, 'height': w + 1, 'method': method})
    payload = {
        'status': 'PASS',
        'claim': 'C_n with cap w (1<=w<=n-2), positive rates: a high recurrent class exists iff n = 1 mod 2w; it is unique',
        'rows': rows,
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    (HERE / 'cycle_classification.json').write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'pairs': len(rows),
                      'high': [(r['n'], r['w'], r.get('states')) for r in rows if r['high_recurrent_classes']]}))


if __name__ == '__main__':
    main()
