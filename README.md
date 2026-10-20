# Bounded swaps in pass-and-swap queues

**Counterexamples to a conjecture of Dorsman and Gardner, and a proof for bipartite swapping graphs with odd caps under positive service rates**

Author: Anonymous. Version of 20 October 2026. This is an unrefereed manuscript with a replayable verification package.

**DOI:** [10.5281/zenodo.23515841](https://doi.org/10.5281/zenodo.23515841) (version 1.0.0-candidate). Source repository: <https://github.com/ipitchford/pass-and-swap-bounded-swaps>. Release page: <https://evidencepress.org/releases/pass-and-swap-bounded-swaps/>. Status and assurance limits: [STATUS.md](STATUS.md), [ASSURANCE.md](ASSURANCE.md). Agents should start at [AI_INDEX.md](AI_INDEX.md).

## In plain English

A *pass-and-swap queue* is a queue in which a job that finishes service may hand its place to a later job it is allowed to swap with. That job then passes the turn further back in the same way, and the last job in the chain leaves. A graph, the *swapping graph*, says which jobs may swap. Dorsman and Gardner (Queueing Systems 107, 2024, [doi:10.1007/s11134-024-09914-1](https://doi.org/10.1007/s11134-024-09914-1)) studied two such queues in a loop and limited each service completion to at most `w` swaps. The order of the jobs defines a direction on every edge of the swapping graph. They proved a simple ("product-form") long-run distribution whenever no directed path in this picture is longer than `w` edges. They conjectured that, if the swapping graph can be coloured with `w+1` colours, the system always ends up in such a situation; that is, every state with a longer path is visited only finitely often.

This package shows that the conjecture is false for every `w >= 2`, and maps out where it fails and where it holds:

* On a cycle of length `n`, if every position is served at a positive rate, a recurrent class with a too-long path exists **exactly when `n = 1 (mod 2w)`**, and it is then the only such class. (For `w >= 2` every cycle satisfies the colouring hypothesis; at `w = 1` odd cycles are outside it.) The smallest examples are the 5-cycle with `w = 2` and the 9-cycle with `w = 2`. The 9-cycle shows that long odd cycles do not rescue the conjecture.
* Two odd cycles joined by a path also give counterexamples. Among all connected graphs with at most 7 vertices, 1355 (graph, cap) pairs violate the conjecture, none of them bipartite and none with `w = 1`.
* In the other direction, the conjecture is **true for bipartite swapping graphs whenever `w` is odd** (in particular `w = 1`), provided every position in the queues receives some service. Exact computation finds no counterexample for even `w` either, on all bipartite graphs with at most 9 vertices; we conjecture that it holds for all `w`.
* On the main counterexample class (the high class of the odd cycle `C_(2w+1)`), equal service rates still give a simple product-form law. With additive per-job rates that vary within a queue, no per-orientation rescaling of the original product weights is stationary on that class, for every `w >= 2`, although the error is numerically small. (The two queues may have different constant rates; low classes keep the ordinary product form for any rates.)

## What is proved and what is certified by computation

"Proved" means a complete written proof in the manuscript; the computations are independent checks. "Certified" means an exact finite computation by the listed program, and that the statement is claimed only for the finite range stated.

| Statement | Status | Where | Check |
|---|---|---|---|
| Binding-move lemma: exact list of orientation changes caused by a capped service; reduction of recurrence to a digraph on orientations | Proved (Lemma 3.3, Cor. 3.4) | Sec. 3 | `census_small_graphs.py` compares it with direct simulation on all graphs with at most 6 vertices |
| Service = (truncated) promotion on linear extensions | Proved (Prop. 3.5) | Sec. 3 | `promotion_check.py` (all graphs with at most 5 vertices and four on 6; binary rotation for w = 2..6) |
| Balanced class on C_{2w+1}: closed for every rate function, recurrent under positive rates, (2w+2)!/(w!)^2 states; conjecture false for all w >= 2 | Proved (Thm 4.1, Cor. 4.2) | Sec. 4 | `balanced_class_check.py` (w = 2..6), `independent_check.py`, `family_check.py` |
| Cycles, under positive rates: high recurrent class iff n = 1 (mod 2w); unique such class (without positive rates: closed set) | Proved (Thm 5.2) | Sec. 5 | `cycle_classification.py` (all n <= 14, all caps; simulation for n <= 9) |
| Two odd cycles joined by a path, cap a+b+L+1 | Proved (Thm 6.1) | Sec. 6.1 | `two_cycle_family.py` (20 members, n <= 13; simulation for n <= 9) |
| Census of all connected graphs, n <= 7 (1355 violating pairs; none bipartite; none at w = 1) | Certified | Sec. 6.2, Table 3 | `census_small_graphs.py` |
| Several jobs per class reduce to one job per class (symmetric lift of the rates) | Proved (Prop. 6.2) | Sec. 6.3 | — |
| Multiplicities: C_5 (all 8 vectors in {1,2,3}^5 with 6-8 jobs, up to symmetry), C_7, C_9 lose their high classes; K_3 with two jobs per class still violates | Certified (finite list) | Table 4 | `multiplicity_check.py` |
| Bipartite graphs, odd w, positive rates: conjecture holds | Proved (Thm 7.1) | Sec. 7 | — |
| Bipartite graphs, every w: no violation for n <= 9 (positive rates) and n <= 7 (head-only service) | Certified; general case is Conjecture 7.2 | Sec. 7 | `bipartite_check.py` (n <= 8), `bipartite_check.py 9` |
| Homogeneous law on B_w; it satisfies the source's partial-balance equations | Proved (Thm 8.2, Prop. 8.3) | Sec. 8.1 | `balanced_class_check.py` (w <= 4), `stationary_details.py` (w = 2, 3) |
| On the balanced class of C_(2w+1), additive rates: orientation-scaled product form holds iff rates are constant within each queue, all w >= 2 | Proved (Thm 8.4) | Sec. 8.2 | `heterogeneous_rates.py` (exact ranks, w = 2, 3, 4), `symbolic_rate_check.py`, `product_form_test.py` |
| Size of the deviation on C_5 (total variation), absorption probabilities, 67/252 | Certified (stated rate values) | Table 5 | `stationary_details.py` |

No formal proof assistant was used.

## Files

* `bounded_swaps_pass_and_swap.pdf`, `.tex`: the manuscript (24 pages). The figures are drawn in TikZ inside the `.tex` file.
* Verifiers (Python 3, standard library only, exact integer/rational arithmetic, deterministic) and their certificates:

| Verifier | Certificate | Run time (one core, approx.) |
|---|---|---|
| `compare_implementations.py` | `transition_comparison.json` | < 1 s |
| `independent_check.py` | `independent_c5.json` | < 1 s |
| `product_form_test.py` | `product_form_rank.json` | < 1 s |
| `symbolic_rate_check.py` | `symbolic_rate_certificate.json` | < 1 s |
| `family_check.py` | `family_head_bijection.json` | < 1 s |
| `balanced_class_check.py` | `balanced_class_check.json` | 8 s |
| `cycle_classification.py` | `cycle_classification.json` | 12 s |
| `two_cycle_family.py` | `two_cycle_family.json` | 22 s |
| `census_small_graphs.py` | `census_small_graphs.json` | 100 s |
| `bipartite_check.py` | `bipartite_check.json` | 160 s |
| `bipartite_check.py 9` | `bipartite_check_n9.json` | about 2 h |
| `multiplicity_check.py` | `multiplicity_check.json` | 1 s |
| `promotion_check.py` | `promotion_check.json` | 2 s |
| `heterogeneous_rates.py` | `heterogeneous_rates.json` | 4 s |
| `stationary_details.py` | `stationary_details.json` | 7 s |

* `verify_all.py`: one-command replay (see below).
* `novelty_report.md`: what is known, what is claimed as new, and the limits of the literature check.
* `research_metadata.json`: machine-readable summary.
* `negative_controls.py`: semantic negative controls. It mutates the capped transition and the claimed cycle rule in temporary copies, and tampers with a certificate and a hashed file, and checks that each change is rejected.
* `SHA256SUMS`: SHA-256 hashes of every other file in the repository.
* Release records: `STATUS.md`, `ASSURANCE.md` (what has and has not been checked), `REVIEW_SUMMARY.md` (substantive objections from internal AI-assisted review and how they were resolved), `REPLAY_RECEIPT.md`, `PROVENANCE.md`, `SOURCES.md`, `CLAIMS.json`, `RESEARCH_GATES.json`, `ENVIRONMENT.txt`, `CITATION.cff`, `.zenodo.json`. The manuscript is the PDF and its LaTeX source.
* Licences: manuscript, prose, certificates and data CC0-1.0 (`LICENSE`, `PUBLIC_DOMAIN.md`); original code MIT (`LICENSE-CODE`); component map in `LICENSES.md`.

## How to replay

```sh
python3 verify_all.py            # every verifier except the 2-hour n = 9 run; about 6 minutes
python3 verify_all.py --long     # including bipartite_check.py 9
python3 verify_all.py --sums     # also check SHA256SUMS
python3 -S -O -B negative_controls.py   # semantic negative controls, about 30 seconds
```

`verify_all.py` copies the scripts into a fresh temporary directory and runs each one there with `python3 -S -O`. `-S` disables site packages, and `-O` shows that no check relies on `assert`: failures raise `RuntimeError`. It then compares every certificate it produced, byte for byte, with the shipped copy. The package directory is never modified. Any single verifier can also be run directly, e.g. `python3 -S -O cycle_classification.py`; it rewrites its certificate next to itself.

To rebuild the PDF (pdfLaTeX with `amsmath`, `amsthm`, `mathtools`, `tikz`, `cleveref`, `aliascnt`, `booktabs`, `tabularx`, `enumitem`, `caption`, `microtype`, `lmodern`, `hyperref`):

```sh
pdflatex bounded_swaps_pass_and_swap.tex
pdflatex bounded_swaps_pass_and_swap.tex
```

## Certificate details

Every new certificate (all except the five inherited ones listed below) contains `input_sha256`, the SHA-256 hash of the verifier source that produced it, and `status: PASS`. It is written with sorted keys and one top-level key per line, so that replays are byte-identical.

* **Model conventions.** The convention is the same in all scripts. Vertices are `0..n-1`. A state is a pair of tuples (upper queue, lower queue), each listed from head to tail. The global word is the upper queue followed by the reversed lower queue. An orientation is stored as a bit mask over the edges `(u,v)`, `u<v`, sorted lexicographically, where a set bit means `u -> v`. Graphs in the census are stored by a canonical upper-triangle code: bit `k` corresponds to the `k`-th pair `(i,j)`, `i<j`, in lexicographic order. They are generated by vertex augmentation with a colour-refinement canonical form, and the counts are compared with OEIS A000088/A001349 (all graphs) and A033995/A005142 (bipartite graphs).
* **Orientation quotient.** Under positive rates all states with one orientation communicate (Lemma 3.2). The census, cycle, family and multiplicity checks therefore compute the digraph of orientation changes. For `n <= 9` (and for every check marked "simulation") they obtain it by simulating the capped transition from every state and every served position. Larger cases use the binding-move lemma, which `census_small_graphs.py` compares with simulation, orientation by orientation, on all 104,590 (graph, cap, orientation) triples with at most 6 vertices. It also compares the closed classes of the full state graph (without any quotient) on 666 (graph, cap) pairs.
* **Inherited certificates.** The five inherited scripts and certificates come from the first version of the work and are unchanged: `compare_implementations.py`, `independent_check.py`, `product_form_test.py`, `symbolic_rate_check.py` and `family_check.py`. They do not record a source hash. `independent_c5.json` also contains an exact stationary law on a nine-orbit quotient of the 180-state class at unit rates. The quotient is exact (lumpable) because the rotations and reflections of `C_5`, together with the exchange of the two queues (which maps a state `(c;d)` to `(d;c)`), map the chain to itself when all rates are one. The script also checks lumpability directly, orbit by orbit. `symbolic_rate_certificate.json` lists every predecessor event of the two states used in Example 8.6, with exact rational-function coefficients in the indeterminate `a`.

## Scope and limits

* One job per vertex unless stated otherwise (Section 6.3 treats several jobs per class), simple swapping graphs, the same cap in both queues.
* "Positive rates" means every position of every queue is served at a positive rate. Closure of the counterexample classes holds for every admissible (order-independent) rate function; recurrence of the whole class, and all transience results, assume positive rates.
* The heterogeneous theorem excludes one representation (original product factors times one constant per orientation). It does not exclude other closed forms. The stationary law exists and is computed exactly for the cases in Table 5.
* The even-cap bipartite case (Conjecture 7.2) and the general head-only-service case remain open.
