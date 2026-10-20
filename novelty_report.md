# Prior work and novelty statement

Version of 20 October 2026. This note records, factually, what the manuscript takes from earlier work, what it claims as new, and how far the literature check went. It contains no assessment of quality or significance.

## The target

Conjecture 1 of J.-P. Dorsman and K. Gardner, *New directions in pass-and-swap queues*, Queueing Systems 107 (2024) 205–256, Section 5.2.2 ([doi:10.1007/s11134-024-09914-1](https://doi.org/10.1007/s11134-024-09914-1), open access under CC BY 4.0). The manuscript quotes the conjecture verbatim together with the source's definitions, and reconciles the conventions in Section 2.3. The reconciliation covers partite versus colourable, the two path-length readings, one job per class, rates, transience and the transition rule.

## Ingredients taken from earlier work (no novelty claimed)

* The pass-and-swap mechanism, invariance of the placement orientation without a cap, and the closed-tandem product form: Comte and Dorsman (2021); Dorsman and Gardner (2024, Lemma 1, Theorems 3 and 7).
* The colouring–longest-path duality (Gallai–Hasse–Roy–Vitaver; Roy 1967; Dorsman and Gardner, Lemma 3).
* Promotion and extended promotion on linear extensions, and the operators tau_i (Schützenberger 1972; Stanley 2009; Ayyer, Klee and Schilling 2014). Reciprocal prefix-sum weights and uniform stationary laws obtained from bijections are known for promotion chains on a fixed poset (Ayyer, Klee and Schilling 2014, Theorems 4.3 and 4.5).
* Cyclic rotation of binary words and its orbit structure (e.g. Reiner, Stanton and White 2004). The scalar sequences that appear are classical: (2w+1)!/(w!)^2 is OEIS A002457, and (2w+2)!/(w!)^2 is an index shift of A002736.
* Source-to-sink reorientations of acyclic orientations (Pretzel 1986; Propp, EJC 2025 / arXiv 2002; Develin, Macauley and Reiner 2016; Defant 2023).

## Claimed as new in the manuscript

1. The binding-move lemma (Lemma 3.3), its reduction of recurrence to a digraph on orientations (Corollary 3.4), and the description of capped service as truncated promotion (Proposition 3.5).
2. The counterexample family on C_{2w+1} for every w >= 2 (Theorem 4.1), together with the observation that closure holds for every admissible rate function. This was in the first version of this work.
3. The classification on cycles (Theorem 5.2): a recurrent class of excessive height exists iff n = 1 (mod 2w).
4. The two-odd-cycle family (Theorem 6.1), the census of connected graphs with at most seven vertices, and the multiplicity reduction and examples (Proposition 6.2, Table 4).
5. The proof of the conjecture for bipartite swapping graphs with odd caps under positive rates (Theorem 7.1), and the even-cap conjecture with its certified range.
6. The homogeneous stationary law on the balanced class and its partial-balance property (Theorem 8.2, Proposition 8.3). The characterisation of when orientation-scaled product weights are stationary (Theorem 8.4) and the quantitative comparison (Table 5) are also new.

## Extent of the literature check

The first version (19 October 2026) checked the full texts of Dorsman and Gardner (2024) and of Comte and Dorsman (2021), and ran public web, arXiv and exact-formula searches. For this revision (20 October 2026) we searched Crossref for works since 2024 matching "pass-and-swap queues", searched arXiv for "pass-and-swap", and queried Semantic Scholar for papers citing Dorsman and Gardner (2024). These searches found no later work on swapping limits and no resolution of the conjecture. They did not include subscription indexes (MathSciNet, zbMATH, Web of Science), non-English sources or correspondence with the authors, so they establish no historical priority. Every reference in the manuscript was checked against its DOI registration (Crossref) or, for arXiv and OEIS, against the primary site.
