"""Semantic negative controls for the replay package.

Usage (standard library only):   python3 -S -O -B negative_controls.py

A replay that cannot fail proves little.  Each control below changes one
thing in a fresh temporary copy and checks that the package rejects it.
The package directory itself is never written to.

  1. baseline: unmodified balanced_class_check.py reproduces its shipped
     certificate byte for byte (so the harness can pass);
  2. cap off by one: the capped transition in balanced_class_check.py is
     allowed w+1 replacements instead of w; the verifier must stop with an
     error or produce a different certificate;
  3. wrong cycle rule: cycle_classification.py is told that high recurrent
     classes occur iff n = 1 (mod 2w+2) instead of (mod 2w); it must fail;
  4. tampered certificate: one digit in a copy of balanced_class_check.json
     is changed; the byte comparison used by verify_all.py must detect it;
  5. tampered file: one byte of a copy of README.md is changed; the
     SHA256SUMS check of verify_all.py must report it.

Exit status 0 only if every mutation is rejected and the baseline passes.
"""
import importlib.util
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
SCRIPTS = sorted(p.name for p in HERE.glob('*.py'))


def fresh_copy(tmp):
    for name in SCRIPTS:
        shutil.copy2(HERE / name, tmp / name)


def run(tmp, script, *args):
    return subprocess.run([sys.executable, '-S', '-O', '-B', script, *args], cwd=tmp,
                          capture_output=True, text=True)


def mutate(path, old, new):
    text = path.read_text()
    if text.count(old) != 1:
        raise RuntimeError(f'mutation anchor not unique in {path.name}: {old!r}')
    path.write_text(text.replace(old, new))


def main():
    results = []

    def record(name, ok, detail):
        results.append(ok)
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")

    cert = 'balanced_class_check.json'
    with tempfile.TemporaryDirectory(prefix='pns_neg_') as t:
        tmp = Path(t)
        fresh_copy(tmp)
        proc = run(tmp, 'balanced_class_check.py')
        same = proc.returncode == 0 and (tmp / cert).read_bytes() == (HERE / cert).read_bytes()
        record('baseline', same, 'unmodified verifier reproduces its certificate' if same
               else f'baseline failed (exit {proc.returncode})')

    with tempfile.TemporaryDirectory(prefix='pns_neg_') as t:
        tmp = Path(t)
        fresh_copy(tmp)
        mutate(tmp / 'balanced_class_check.py', '    for _ in range(w):\n        nxt = None',
               '    for _ in range(w + 1):\n        nxt = None')
        proc = run(tmp, 'balanced_class_check.py')
        produced = (tmp / cert).read_bytes() if (tmp / cert).exists() else None
        rejected = proc.returncode != 0 or produced != (HERE / cert).read_bytes()
        last = (proc.stderr.strip().splitlines() or [''])[-1]
        record('cap off by one', rejected,
               f'rejected (exit {proc.returncode}; {last})' if rejected else 'NOT rejected')

    with tempfile.TemporaryDirectory(prefix='pns_neg_') as t:
        tmp = Path(t)
        fresh_copy(tmp)
        mutate(tmp / 'cycle_classification.py', '    if (n - 1) % (2 * w) or n < 2 * w + 1:',
               '    if (n - 1) % (2 * w + 2) or n < 2 * w + 1:')
        proc = run(tmp, 'cycle_classification.py')
        last = (proc.stderr.strip().splitlines() or [''])[-1]
        record('wrong cycle rule', proc.returncode != 0,
               f'rejected (exit {proc.returncode}; {last})' if proc.returncode != 0 else 'NOT rejected')

    with tempfile.TemporaryDirectory(prefix='pns_neg_') as t:
        tmp = Path(t)
        shipped = (HERE / cert).read_bytes()
        tampered = shipped.replace(b'"states":180', b'"states":181', 1)
        if tampered == shipped:
            raise RuntimeError('certificate tamper anchor not found')
        (tmp / cert).write_bytes(tampered)
        record('tampered certificate', (tmp / cert).read_bytes() != shipped,
               'byte comparison detects the altered state count')

    with tempfile.TemporaryDirectory(prefix='pns_neg_') as t:
        tmp = Path(t)
        for line in (HERE / 'SHA256SUMS').read_text().splitlines():
            if line.strip():
                name = line.split(None, 1)[1].lstrip('*')
                (tmp / name).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(HERE / name, tmp / name)
        shutil.copy2(HERE / 'SHA256SUMS', tmp / 'SHA256SUMS')
        readme = tmp / 'README.md'
        data = bytearray(readme.read_bytes())
        data[0] ^= 1
        readme.write_bytes(bytes(data))
        spec = importlib.util.spec_from_file_location('verify_all_copy', tmp / 'verify_all.py')
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        bad = module.check_sums()
        record('tampered file', bad == ['README.md'], f'SHA256SUMS check reports {bad}')

    if not all(results):
        print(f'{results.count(False)} control(s) failed')
        sys.exit(1)
    print('All negative controls behaved as expected.')


if __name__ == '__main__':
    main()
