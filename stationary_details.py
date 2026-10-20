"""Quantitative details on the five-cycle class B_2 and partial balance.

Certificate: stationary_details.json.  Standard library only; exact
rational arithmetic; deterministic.  Run:  python3 -S -O stationary_details.py

  (1) For rates (1,1,1,1,a) at both queues, a in {1/2, 2, 10}: the exact
      stationary law pi on B_2 (180 states) is computed by exact Gaussian
      elimination and compared with
        - q_r normalised on B_2, and
        - the orientation-matched weights pi(Sigma_A) q_r / q_r(Sigma_A).
      Total-variation distances, the spread of pi/q_r within orientations,
      the mean upper-queue length and the upper-queue throughput are recorded.
  (2) Unit rates, full 720-state chain on C_5 with cap 2: exact probability of
      absorption into B_2 from every transient state, and the long-run mass of
      B_2 from the uniform initial distribution.
  (3) Homogeneous rates mu0 = 2, nu0 = 3, caps w = 2, 3: the law of the paper
      satisfies the partial-balance equations (23)-(24) of Dorsman--Gardner
      (2024) at every state of B_w (outflow by upper services = inflow by
      lower services, and conversely).
"""
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations, permutations
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


def okey(word, n):
    pos = {v: i for i, v in enumerate(word)}
    return tuple(int(pos[v] < pos[(v + 1) % n]) for v in range(n))


def phi(seq, rate):
    out = F(1)
    tot = F(0)
    for v in seq:
        tot += rate[v]
        out /= tot
    return out


def events(state, adj, w, mu, nu):
    c, d = state
    out = []
    for p in range(len(c)):
        nc, dep = capped_transfer(c, p, adj, w)
        out.append(((nc, d + (dep,)), mu[c[p]], 'up'))
    for p in range(len(d)):
        nd, dep = capped_transfer(d, p, adj, w)
        out.append(((c + (dep,), nd), nu[d[p]], 'low'))
    return out


def solve(A, b):
    """Exact Gauss-Jordan for a square nonsingular system."""
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for col in range(n):
        piv = next(i for i in range(col, n) if M[i][col])
        M[col], M[piv] = M[piv], M[col]
        f = M[col][col]
        M[col] = [x / f for x in M[col]]
        for i in range(n):
            if i != col and M[i][col]:
                g = M[i][col]
                M[i] = [x - g * y for x, y in zip(M[i], M[col])]
    return [M[i][n] for i in range(n)]


def stationary(states, adj, w, mu, nu):
    idx = {s: i for i, s in enumerate(states)}
    N = len(states)
    Q = [[F(0)] * N for _ in range(N)]
    for i, s in enumerate(states):
        for t, rate, _ in events(s, adj, w, mu, nu):
            j = idx[t]
            Q[i][j] += rate
            Q[i][i] -= rate
    # pi Q = 0, sum pi = 1: transpose, replace last equation by normalisation
    A = [[Q[i][j] for i in range(N)] for j in range(N)]
    A[-1] = [F(1)] * N
    b = [F(0)] * (N - 1) + [F(1)]
    pi = solve(A, b)
    for j in range(N):
        check(sum(pi[i] * Q[i][j] for i in range(N)) == 0, 'not stationary')
    return pi


def part1():
    w, n = 2, 5
    adj = [(1 << ((v + 1) % n)) | (1 << ((v - 1) % n)) for v in range(n)]
    states = [(W[:k], tuple(reversed(W[k:]))) for W in balanced_words(w) for k in range(n + 1)]
    rows = []
    for a in (F(1, 2), F(2), F(10)):
        r = [F(1)] * 4 + [a]
        pi = stationary(states, adj, w, r, r)
        q = [phi(c, r) * phi(d, r) for c, d in states]
        Z = sum(q)
        qn = [x / Z for x in q]
        ok = [okey(c + tuple(reversed(d)), n) for c, d in states]
        massp = {}
        massq = {}
        for i, o in enumerate(ok):
            massp[o] = massp.get(o, F(0)) + pi[i]
            massq[o] = massq.get(o, F(0)) + q[i]
        matched = [massp[ok[i]] * q[i] / massq[ok[i]] for i in range(len(states))]
        tv_plain = sum(abs(x - y) for x, y in zip(pi, qn)) / 2
        tv_matched = sum(abs(x - y) for x, y in zip(pi, matched)) / 2
        spread = []
        for o in sorted(massp):
            ratios = [pi[i] / q[i] for i in range(len(states)) if ok[i] == o]
            spread.append(max(ratios) / min(ratios))
        mean_c = sum(pi[i] * len(states[i][0]) for i in range(len(states)))
        mean_c_m = sum(matched[i] * len(states[i][0]) for i in range(len(states)))
        thr = sum(pi[i] * sum((r[v] for v in states[i][0]), F(0)) for i in range(len(states)))
        thr_m = sum(matched[i] * sum((r[v] for v in states[i][0]), F(0)) for i in range(len(states)))
        rows.append({'a': str(a), 'tv_to_normalised_q': str(tv_plain), 'tv_to_normalised_q_float': float(tv_plain),
                     'tv_to_orientation_matched_q': str(tv_matched),
                     'tv_to_orientation_matched_q_float': float(tv_matched),
                     'max_over_orientations_of_max_over_min_pi_over_q': float(max(spread)),
                     'mean_upper_length_pi': str(mean_c), 'mean_upper_length_matched': str(mean_c_m),
                     'upper_throughput_pi': str(thr), 'upper_throughput_matched': str(thr_m),
                     'orientation_masses_pi': sorted(str(x) for x in massp.values())})
    return rows


def part2():
    w, n = 2, 5
    adj = [(1 << ((v + 1) % n)) | (1 << ((v - 1) % n)) for v in range(n)]
    one = [F(1)] * n
    states = [(p[:k], tuple(reversed(p[k:]))) for p in permutations(range(n)) for k in range(n + 1)]
    idx = {s: i for i, s in enumerate(states)}
    g = {i: {idx[t] for t, _, _ in events(s, adj, w, one, one)} for i, s in enumerate(states)}
    bal = {(W[:k], tuple(reversed(W[k:]))) for W in balanced_words(w) for k in range(n + 1)}
    bal_i = {idx[s] for s in bal}
    closed = []
    for comp in sccs(g):
        cs = set(comp)
        if all(t in cs for s in comp for t in g[s]):
            closed.append(cs)
    check(sorted(len(c) for c in closed) == [48] * 10 + [180], 'closed classes')
    check(any(c == bal_i for c in closed), 'B_2 is not a closed class')
    rec = set().union(*closed)
    trans = sorted(set(range(len(states))) - rec)
    check(len(trans) == 60, 'transient count')
    ti = {s: k for k, s in enumerate(trans)}
    # h = P h on transient states, h = 1 on B_2, 0 on other closed classes
    A = [[F(0)] * len(trans) for _ in trans]
    b = [F(0)] * len(trans)
    for k, s in enumerate(trans):
        A[k][k] += 1
        evs = events(states[s], adj, w, one, one)
        tot = sum(r for _, r, _ in evs)
        for t, r, _ in evs:
            j = idx[t]
            if j in ti:
                A[k][ti[j]] -= r / tot
            elif j in bal_i:
                b[k] += r / tot
    h = solve(A, b)
    values = {}
    for x in h:
        values[str(x)] = values.get(str(x), 0) + 1
    mass = (F(len(bal_i)) + sum(h)) / len(states)
    check(mass == F(67, 252), 'uniform-start mass of B_2 should be 67/252')
    return {'absorption_probability_counts': values, 'uniform_start_long_run_mass_of_B2': str(mass)}


def part3():
    out = []
    for w in (2, 3):
        n = 2 * w + 1
        adj = [(1 << ((v + 1) % n)) | (1 << ((v - 1) % n)) for v in range(n)]
        mu0, nu0 = F(2), F(3)
        mu = [mu0] * n
        nu = [nu0] * n
        states = [(W[:k], tuple(reversed(W[k:]))) for W in balanced_words(w) for k in range(n + 1)]
        sset = set(states)
        M = len(states) // (n + 1)
        p = nu0 / (mu0 + nu0)
        pi = {s: F(comb(n, len(s[0]))) * p ** len(s[0]) * (1 - p) ** (n - len(s[0])) / M for s in states}
        check(sum(pi.values()) == 1, 'normalisation')
        inflow_from_lower = {s: F(0) for s in states}
        inflow_from_upper = {s: F(0) for s in states}
        for s in states:
            for t, r, kind in events(s, adj, w, mu, nu):
                check(t in sset, 'not closed')
                if kind == 'low':
                    inflow_from_lower[t] += pi[s] * r
                else:
                    inflow_from_upper[t] += pi[s] * r
        for s in states:
            c, d = s
            if c:
                check(pi[s] * mu0 * len(c) == inflow_from_lower[s], 'partial balance (23) fails')
            if d:
                check(pi[s] * nu0 * len(d) == inflow_from_upper[s], 'partial balance (24) fails')
        out.append({'w': w, 'states': len(states), 'mu0': '2', 'nu0': '3', 'partial_balance_23_24': True})
    return out


def main():
    payload = {
        'status': 'PASS',
        'heterogeneous_c5_size_of_deviation': part1(),
        'absorption_unit_rates_c5_w2': part2(),
        'homogeneous_partial_balance': part3(),
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    (HERE / 'stationary_details.json').write_text(dump(payload))
    print(json.dumps({'status': 'PASS'}))


if __name__ == '__main__':
    main()
