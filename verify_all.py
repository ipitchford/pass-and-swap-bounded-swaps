"""Replay every verifier of this package and compare certificates byte for byte.

Usage (from any directory; Python 3.9 or later, standard library only):

    python3 verify_all.py            # all verifiers except the long n = 9 run (about 6 minutes)
    python3 verify_all.py --long     # also bipartite_check.py 9 (about 2 CPU hours)
    python3 verify_all.py --sums     # additionally check SHA256SUMS

Each verifier is copied, together with the other scripts, into a fresh
temporary directory and run there as  python3 -S -O <script> [args].  The
certificate it writes is compared with the shipped certificate of the same
name.  The package directory itself is never written to.  The exit status is
0 only if every run succeeds and every certificate is identical.
"""
import hashlib
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

RUNS = [
    # (script, arguments, certificates written, long?)
    ('compare_implementations.py', [], ['transition_comparison.json'], False),
    ('independent_check.py', [], ['independent_c5.json'], False),
    ('product_form_test.py', [], ['product_form_rank.json'], False),
    ('symbolic_rate_check.py', [], ['symbolic_rate_certificate.json'], False),
    ('family_check.py', [], ['family_head_bijection.json'], False),
    ('balanced_class_check.py', [], ['balanced_class_check.json'], False),
    ('cycle_classification.py', [], ['cycle_classification.json'], False),
    ('two_cycle_family.py', [], ['two_cycle_family.json'], False),
    ('census_small_graphs.py', [], ['census_small_graphs.json'], False),
    ('bipartite_check.py', [], ['bipartite_check.json'], False),
    ('bipartite_check.py', ['9'], ['bipartite_check_n9.json'], True),
    ('multiplicity_check.py', [], ['multiplicity_check.json'], False),
    ('promotion_check.py', [], ['promotion_check.json'], False),
    ('heterogeneous_rates.py', [], ['heterogeneous_rates.json'], False),
    ('stationary_details.py', [], ['stationary_details.json'], False),
]


def check_sums():
    sums = HERE / 'SHA256SUMS'
    bad = []
    for line in sums.read_text().splitlines():
        if not line.strip():
            continue
        digest, name = line.split(None, 1)
        name = name.lstrip('*')
        actual = hashlib.sha256((HERE / name).read_bytes()).hexdigest()
        if actual != digest:
            bad.append(name)
    return bad


def main():
    long_runs = '--long' in sys.argv
    failures = []
    scripts = sorted({r[0] for r in RUNS})
    for script, args, certs, is_long in RUNS:
        if is_long and not long_runs:
            print(f'SKIP  {script} {" ".join(args)}  (use --long)')
            continue
        with tempfile.TemporaryDirectory(prefix='pns_replay_') as tmp:
            tmp = Path(tmp)
            for s in scripts:
                shutil.copy2(HERE / s, tmp / s)
            start = time.time()
            proc = subprocess.run([sys.executable, '-S', '-O', script, *args], cwd=tmp,
                                  capture_output=True, text=True)
            elapsed = time.time() - start
            if proc.returncode != 0:
                failures.append(f'{script}: exit status {proc.returncode}')
                print(f'FAIL  {script} {" ".join(args)}  exit {proc.returncode}\n{proc.stderr[-2000:]}')
                continue
            for cert in certs:
                new = (tmp / cert).read_bytes() if (tmp / cert).exists() else None
                old = (HERE / cert).read_bytes()
                if new != old:
                    failures.append(f'{script}: {cert} differs')
                    print(f'FAIL  {script} {" ".join(args)}  {cert} differs from the shipped copy')
                else:
                    print(f'PASS  {script} {" ".join(args)}  {cert} identical  ({elapsed:.0f} s)')
    if '--sums' in sys.argv:
        bad = check_sums()
        if bad:
            failures.append('SHA256SUMS mismatch: ' + ', '.join(bad))
            print('FAIL  SHA256SUMS:', ', '.join(bad))
        else:
            print('PASS  SHA256SUMS')
    if failures:
        print(f'{len(failures)} failure(s)')
        sys.exit(1)
    print('All replays identical.')


if __name__ == '__main__':
    main()
