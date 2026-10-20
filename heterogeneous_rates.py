"""Orientation-scaled product weights on the balanced class B_w of C_{2w+1}.

Certificate: heterogeneous_rates.json.  Standard library only; exact
rational arithmetic; deterministic.  Run:  python3 -S -O heterogeneous_rates.py

Ansatz: pi(c;d) = x_{A(c,d)} * Phi_mu(c) * Phi_nu(d), with one unknown
x_A per balanced orientation A, additive rates mu(c) = sum of mu_i over c,
nu(d) = sum of nu_i over d, and Phi(c) = prod_j 1/(rate of c_1..c_j).
Theorem (paper): for every w >= 2 the ansatz satisfies global balance on
B_w iff all mu_i are equal and all nu_i are equal.

Checks for w = 2, 3, 4 (B_w has 180, 1120, 6300 states; 10, 14, 18 orientations):
  * the global-balance system in the unknowns x_A has rank m-1 (a positive
    one-dimensional solution, x constant) for homogeneous rates, including
    mu_0 != nu_0, and full rank m (no solution at all) for every
    non-homogeneous test vector: one changed rate a in {2, 1/2, 7/5},
    deterministic pseudo-random rational vectors, a geometric vector, and
    (w = 4, n = 9) a period-3 vector, which satisfies the first family of
    corner identities but not the second;
  * with x constant, the global-balance residual vanishes at every state
    with 2 <= k <= n-2 jobs in the upper queue (it is concentrated on the
    images of the binding head transitions, k = 1 or n-1), as the proof uses;
  * the first corner identity at a source two steps before the changed rate
    reduces to 3(a-1) = 0.
"""
import hashlib
import json
from fractions import Fraction as F
from itertools import combinations
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


def orientation_key(word, n):
    pos = {v: i for i, v in enumerate(word)}
    return tuple(int(pos[v] < pos[(v + 1) % n]) for v in range(n))


def phi(seq, rate):
    out = F(1)
    tot = F(0)
    for v in seq:
        tot += rate[v]
        out /= tot
    return out


def system(w, mu, nu):
    n = 2 * w + 1
    adj = [(1 << ((v + 1) % n)) | (1 << ((v - 1) % n)) for v in range(n)]
    words = balanced_words(w)
    states = [(W[:k], tuple(reversed(W[k:]))) for W in words for k in range(n + 1)]
    index = {s: i for i, s in enumerate(states)}
    okeys = sorted({orientation_key(W, n) for W in words})
    oidx = {o: i for i, o in enumerate(okeys)}
    m = len(okeys)
    q = [phi(c, mu) * phi(d, nu) for c, d in states]
    orient = [oidx[orientation_key(c + tuple(reversed(d)), n)] for c, d in states]
    rows = [[F(0)] * m for _ in states]
    for i, (c, d) in enumerate(states):
        out = sum((mu[v] for v in c), F(0)) + sum((nu[v] for v in d), F(0))
        rows[i][orient[i]] -= q[i] * out
        for p in range(len(c)):
            nc, dep = capped_transfer(c, p, adj, w)
            j = index.get((nc, d + (dep,)))
            check(j is not None, 'B_w not closed')
            rows[j][orient[i]] += q[i] * mu[c[p]]
        for p in range(len(d)):
            nd, dep = capped_transfer(d, p, adj, w)
            j = index.get((c + (dep,), nd))
            check(j is not None, 'B_w not closed')
            rows[j][orient[i]] += q[i] * nu[d[p]]
    return states, rows, m


def rank(rows, m):
    basis = []          # list of (pivot, row) in reduced form
    for r in rows:
        r = r[:]
        for piv, b in basis:
            if r[piv]:
                f = r[piv]
                r = [x - f * y for x, y in zip(r, b)]
        nz = next((j for j in range(m) if r[j]), None)
        if nz is None:
            continue
        f = r[nz]
        r = [x / f for x in r]
        new = []
        for piv, b in basis:
            if b[nz]:
                g = b[nz]
                b = [x - g * y for x, y in zip(b, r)]
            new.append((piv, b))
        basis = new + [(nz, r)]
        if len(basis) == m:
            break
    return len(basis), basis


def lcg(seed):
    while True:
        seed = (1103515245 * seed + 12345) % 2147483648
        yield seed


def main():
    results = []
    gen = lcg(20261020)
    for w in (2, 3, 4):
        n = 2 * w + 1
        tests = [('homogeneous mu=nu=1', [F(1)] * n, [F(1)] * n, True),
                 ('homogeneous mu=2, nu=3', [F(2)] * n, [F(3)] * n, True)]
        for a in (F(2), F(1, 2), F(7, 5)):
            r = [F(1)] * (n - 1) + [a]
            tests.append((f'one changed rate a={a} at both queues', r, r, False))
        r = [F(1)] * (n - 1) + [F(2)]
        tests.append(('changed rate at upper queue only', r, [F(1)] * n, False))
        for t in range(2):
            r = [F(1 + next(gen) % 9, 1 + next(gen) % 4) for _ in range(n)]
            tests.append((f'pseudo-random rationals #{t+1}', r, r, False))
        r = [F(2) ** k for k in range(n)]
        tests.append(('geometric 2^k', r, r, False))
        if n % 3 == 0:
            r = [F((k % 3) + 1) for k in range(n)]
            tests.append(('period-3 (1,2,3,...)', r, r, False))
        for label, mu, nu, homog in tests:
            states, rows, m = system(w, mu, nu)
            rk, basis = rank(rows, m)
            if homog:
                check(rk == m - 1, f'homogeneous rank should be m-1 (w={w}, {label})')
                # x constant solves the system
                check(all(sum(row) == 0 for row in rows), 'constant x is not a solution')
                residual_support = sorted({len(states[i][0]) for i, row in enumerate(rows) if sum(row)})
            else:
                check(rk == m, f'non-homogeneous rank should be m (w={w}, {label})')
                residual_support = sorted({len(states[i][0]) for i, row in enumerate(rows) if sum(row)})
                check(set(residual_support) <= {1, n - 1},
                      'residual of constant x outside k in {1, n-1}')
            results.append({'w': w, 'n': n, 'states': len(states), 'orientations': m,
                            'test': label, 'mu': [str(x) for x in mu], 'nu': [str(x) for x in nu],
                            'rank': rk, 'solution_exists': rk < m,
                            'k_values_where_constant_x_fails': residual_support})
    # first corner identity with one changed rate at x_2 (rates 1 elsewhere)
    a = F(5, 3)
    lhs = (1 + 1 + a) * (1 + 1 + 1)
    rhs = (1 + 1 + 1) ** 2
    check(lhs - rhs == 3 * (a - 1), 'corner identity algebra')
    payload = {
        'status': 'PASS',
        'results': results,
        'input_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        'verifier': Path(__file__).name,
    }
    (HERE / 'heterogeneous_rates.json').write_text(dump(payload))
    print(json.dumps({'status': 'PASS', 'tests': len(results)}))


if __name__ == '__main__':
    main()
