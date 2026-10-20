# Review summary (internal, AI-assisted)

This is a public summary of the internal checks applied before release. It keeps the substantive objections and how each was resolved. It is **not** peer review: all reviewers were AI systems working inside the Evidence Press workflow, and none of this is independent reproduction by an unaffiliated party. The original review reports are kept privately; this summary contains no scores or ratings.

## What was done

1. **First round (19-20 October 2026).** A multi-perspective AI-assisted review of the first version: an editor's report, three specialist perspectives (Markov-chain combinatorics and exact verification; product-form and pass-and-swap queues; algebraic combinatorics and applications) and a devil's-advocate challenge to the core argument. Each reviewer wrote its own code and reproduced the C_5 census (720 states: ten closed 48-state classes, one closed 180-state class, 60 transient states) and the balanced-class sizes for larger caps. No reviewer found a mathematical error in the counterexample.
2. **Revision.** The manuscript was substantially rewritten and extended in response (see below).
3. **Verification re-review.** A re-review that wrote a separate engine from Definition 2.2 without using the package scripts, read the new proofs line by line and recomputed every new numerical claim (Tables 2-5, the C_9 class with 131,040 states, the census rows, the multiplicity examples, partial balance and the nullity computations behind Theorem 8.4). It found no mathematical error and listed wording and scope corrections.
4. **Adversarial review by a different model.** One AI-assisted adversarial review (an OpenAI reasoning model, 20 October 2026) was asked for the strongest counterexample or gap in each principal result. It did not refute any principal result and raised scope points, all addressed below.
5. **Clean replay.** All 15 verifier runs, including the 9-vertex bipartite census, were replayed from a separate copy and reproduced every certificate byte for byte.

## Substantive objections and their resolution

| Objection | Resolution in this version |
|---|---|
| The conjecture was paraphrased, not quoted, and the source's conventions (partite versus colourable, path length counted in vertices or edges, one job per class, rates) were not reconciled. | Section 2.3 quotes Conjecture 1 and the source's definitions verbatim and reconciles each convention. The counterexamples violate both path-length readings; one job per class follows from the source's Remark 9 and Proposition 6.2. |
| The paper did not explain why the conjecture mattered or give an informal model. | New Section 1: order-independent queues and product form, pass-and-swap queues, the redundancy-scheduling motivation for caps, a worked C_5 example with a figure, a contributions list, a roadmap and a notation table. |
| The heterogeneous-rates result was proved only for C_5, w = 2, one rate pattern and one ansatz, but was described as if product form failed in general. | Theorem 8.4 now covers every w >= 2 and is stated as the exclusion of one representation (original product weights times one constant per orientation), for additive per-job rates. Table 5 reports the size of the deviation (total variation below 2e-3 after per-orientation rescaling on C_5) and the dependence on the initial state (probability 67/252 of reaching the high class from a uniform start). |
| One counterexample family is not a characterisation; reviewers' own computations showed failures beyond odd cycles of length 2w+1. | Added the exact cycle classification (Theorem 5.2: n = 1 mod 2w, which corrected the reviewers' initial guess of n = 2w+1 only), the two-odd-cycle family (Theorem 6.1), the census of all connected graphs with at most 7 vertices, the multiplicity results, and the positive bipartite theorem for odd caps (Theorem 7.1). |
| Literature too thin. | References expanded from 4 to 24 and checked against Crossref or primary sites. |
| Implicit proof steps (carrier path, initiating state, specialisation of the Comte-Dorsman stationary theorem, a quotient used in one check). | Made explicit in Lemma 3.3(b), Remark 4.4 and Theorem 8.2; the quotient is explained in the README. |
| (Re-review) Table 4, second row, misstated the range of multiplicity vectors. | Corrected to "every vector in {1,2,3}^5 with 6 to 8 jobs, up to rotation and reflection (8 vectors)". |
| (Re-review and adversarial review) The abstract and contributions omitted the positive-rates hypothesis for the cycle classification, said the conjecture "holds exactly when n is not 1 mod 2w" without excluding w = 1 (where odd cycles lie outside the hypothesis), and omitted that Theorem 8.4 concerns additive per-job rates. | Qualifiers added in the abstract, Contribution C3, Section 5 and the Discussion. |
| (Re-review) Two sentences over-generalised: one implied failure on every non-bipartite graph; one stated a census result as a general claim. | Both reworded to their proved or computed scope. |
| (Re-review) Proposition 3.5(b) was ambiguous about which job moves at the cap. | Restated: the larger element of the w-th comparable pair (the w-th displaced job) moves to the end. |
| (Adversarial review) Proposition 6.2 needs the symmetric lift of the rate functions as an explicit hypothesis. | The hypothesis is stated in the proposition. |
| (Adversarial review) The lower-queue case of Proposition 3.5 was too terse. | Part (c) now states the queue-exchange involution explicitly. |
| (Adversarial review) The product-form consequence of the conjecture must be read classwise, because the chain is not irreducible. | Section 1.2 states the classwise reading. |
| (Adversarial review) Finite certificates cannot be audited from the manuscript alone. | The verifiers, certificates, one-command replay and hashes are archived under the DOI; the manuscript describes what hashes and replay do and do not establish (Section 10). |
| (Re-review) Agreement of a mean queue length in Table 5 is forced by symmetry and is not evidence of accuracy. | The text now says so. |

## Remaining limitations

* Conjecture 7.2 (bipartite swapping graphs, even caps) is open; it is certified for bipartite graphs with at most 9 vertices.
* Head-only (general order-independent) service is certified only for bipartite graphs with at most 7 vertices.
* There is no general characterisation of the pairs (G, w) with a recurrent class of excessive height.
* The heterogeneous stationary law on the balanced class has no known closed form.
* No proof assistant was used; no unaffiliated party has reproduced the results; the novelty search was bounded.
