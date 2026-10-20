# Replay receipt

All verifier scripts (`*.py` other than `negative_controls.py`) and all JSON certificates (other than `research_metadata.json`) in this release are byte-identical to the reviewed candidate package of 20 October 2026. Publication changed only the DOI text in the manuscript, README and metadata, the rebuilt PDF, the added release records and `negative_controls.py`.

## Runs

| Run | Command | Result |
|---|---|---|
| Producer replay, candidate package, 20 October 2026 | `python3 verify_all.py --long --sums` | All 15 runs byte-identical, including `bipartite_check.py 9`; SHA256SUMS passed |
| Clean replay from a separate copy of the package, 20 October 2026 | `python3 verify_all.py --long --sums` | All 15 runs byte-identical, including `bipartite_check.py 9` (5333 s); SHA256SUMS passed |
| Release replay, this repository, 20 October 2026, CPython 3.14.7, macOS arm64 | `python3 -S -O -B verify_all.py --sums` | 14 runs byte-identical (the 9-vertex run was not repeated because no verifier or certificate changed); SHA256SUMS passed |
| Negative controls, same host | `python3 -S -O -B negative_controls.py` | Baseline reproduces its certificate; off-by-one cap rejected (`RuntimeError: B_w is not strongly connected`); wrong cycle rule rejected (`RuntimeError: unexpected high class n=3 w=1`); tampered certificate and tampered file detected |
| Continuous integration | `.github/workflows/replay.yml` (CPython 3.12 and 3.14) | See the GitHub Actions run for the tagged commit |

Release replay timings (one core): balanced_class_check 6 s, cycle_classification 8 s, two_cycle_family 16 s, census_small_graphs 81 s, bipartite_check 137 s, all others under 8 s.

## What this establishes

Byte-identical replay shows that the shipped certificates are exactly what the shipped programs produce, under `-O` (so no check relies on `assert`). It does not show that the programs implement Definition 2.2 correctly; that rests on the cross-checks described in the README and Section 10 of the manuscript, on the negative controls, and on an internal re-review that wrote its own engine. None of these runs is an independent reproduction by an unaffiliated party.
