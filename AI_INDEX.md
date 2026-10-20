# Agent navigation

## Identity

- Release slug: `pass-and-swap-bounded-swaps`; version `1.0.0-candidate`; DOI 10.5281/zenodo.23515841
- Scholarly creator: Anonymous. Status: unrefereed candidate ([STATUS.md](STATUS.md))
- Manuscript: [PDF](bounded_swaps_pass_and_swap.pdf), [LaTeX source](bounded_swaps_pass_and_swap.tex)
- Machine-readable claims: [CLAIMS.json](CLAIMS.json); summary metadata: [research_metadata.json](research_metadata.json)
- Assurance limits: [ASSURANCE.md](ASSURANCE.md); internal review summary: [REVIEW_SUMMARY.md](REVIEW_SUMMARY.md)
- Provenance: [PROVENANCE.md](PROVENANCE.md); sources: [SOURCES.md](SOURCES.md), [novelty_report.md](novelty_report.md), [RESEARCH_GATES.json](RESEARCH_GATES.json)
- Citation: [CITATION.cff](CITATION.cff), [Zenodo metadata](.zenodo.json); licences: [LICENSES.md](LICENSES.md)

## Scientific map

The object is a closed tandem of two pass-and-swap queues whose service completions may trigger at most `w` swaps (Dorsman and Gardner, Queueing Systems 107, 2024). The relative order of jobs defines an acyclic orientation of the swapping graph; an orientation is *high* if it has a directed path with more than `w` edges.

- Conjecture 1 of Dorsman and Gardner (if the swapping graph is (w+1)-colourable, high states are transient) is **false for every w >= 2**: Theorem 4.1 (balanced class of the odd cycle C_{2w+1}).
- Cycles, under positive rates: a high recurrent class exists iff `n = 1 (mod 2w)`, and it is unique (Theorem 5.2).
- Non-cycle counterexamples: two odd cycles joined by a path (Theorem 6.1); census of all connected graphs with at most 7 vertices (1355 violating pairs, none bipartite, none at w = 1).
- Positive result: bipartite swapping graphs, odd `w`, positive rates: the conjecture holds (Theorem 7.1). Even `w` is Conjecture 7.2 (certified for at most 9 vertices).
- Stationary laws on the counterexample class (Section 8): explicit law for homogeneous rates; Theorem 8.4 excludes orientation-scaled product weights for additive rates that vary within a queue.

Do not read "certified" results as proofs beyond their finite range, Theorem 8.4 as a failure of every closed form, or the cycle result at `w = 1` as a counterexample (odd cycles are outside the hypothesis there).

## Replay

From the repository root (standard-library Python 3.9+):

```sh
python3 -S -O -B verify_all.py --sums      # all runs except the 9-vertex census; about 6 minutes
python3 -S -O -B verify_all.py --long      # adds bipartite_check.py 9; about 1.5-2 CPU hours
python3 -S -O -B negative_controls.py      # semantic negative controls
```

[verify_all.py](verify_all.py) reruns each verifier in a temporary directory and compares certificates byte for byte. [negative_controls.py](negative_controls.py) shows the replay rejects a wrong cap, a wrong cycle rule and tampered files. Results: [REPLAY_RECEIPT.md](REPLAY_RECEIPT.md). Environment: [ENVIRONMENT.txt](ENVIRONMENT.txt). Hashes: [SHA256SUMS](SHA256SUMS).

## Useful next research question

Prove Conjecture 7.2 (bipartite swapping graphs with even caps), or characterise all pairs (G, w) with a recurrent class of excessive height.
