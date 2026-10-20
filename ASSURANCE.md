# Assurance

This file states what has and has not been checked. Each line is a separate kind of assurance; none implies another.

| Dimension | State | Evidence |
|---|---|---|
| Availability | Public GitHub release and Zenodo DOI 10.5281/zenodo.23515841 | Release assets and SHA-256 sums |
| Internal replay | Passed | All 15 verifier runs reproduce their shipped certificates byte for byte (`verify_all.py --long --sums`), including the 9-vertex bipartite census (about 1.5 CPU hours). A second clean replay from a separate copy of the package also passed. See [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). |
| Semantic negative controls | Passed | `negative_controls.py`: an off-by-one cap and a wrong cycle rule are rejected; a tampered certificate and a tampered file are detected. |
| Implementation diversity | Partial | The five inherited checkers (`compare_implementations.py`, `independent_check.py`, ...) were written separately from the main verifiers; the census compares three descriptions (full state graph, simulated orientation digraph, Lemma 3.3) on all graphs with at most 6 vertices. This is diversity inside one producer workflow, not independent reimplementation. |
| Internal AI-assisted review | Done | Multi-perspective AI-assisted review, a verification re-review with its own code, and one AI-assisted adversarial review by a different model; summary in [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md). Not peer review. |
| Independent rerun or reimplementation by an unaffiliated party | Not established | |
| Specialist (human) review | Not established | |
| Editorial peer review | Not established | |
| Formal verification | Not attempted | No proof assistant was used. |
| Environment reproducibility | Standard-library Python 3 only, exact integer and rational arithmetic, deterministic output | [ENVIRONMENT.txt](ENVIRONMENT.txt) |
| Novelty | Bounded search only | [novelty_report.md](novelty_report.md), [RESEARCH_GATES.json](RESEARCH_GATES.json) |

Hashes and byte-identical replay establish integrity and reproducibility of the computations. They do not establish that the programs implement Definition 2.2 of the manuscript correctly; the cross-checks above are the evidence for that. The general theorems rest on the written proofs, not on the computations.
