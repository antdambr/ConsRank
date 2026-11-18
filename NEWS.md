# ConsRank 3.1.2

## Improvements

* `bolt()` (and `consrank(algorithm = "BOLT")`) with ties allowed
  (`full = FALSE`): a first phase runs the local search without ties from each
  starting solution and adds the resulting full rankings to the starting
  solutions. The result is never worse than before (all the previous starting
  solutions are still used) and never worse than `bolt(full = TRUE)` (a full
  ranking is also a bucket order). On data with nearly random judges the
  optimum has few ties, while the previous starting solutions could lead to
  local optima with too many ties: on 707 simulated data sets with certified
  optimum (M = 50 to 200 objects) the optimum is now found in 426 cases instead
  of 391. `bolt(full = TRUE)` and `mcbo()` are unchanged.

# ConsRank 3.1.1

## Bug fixes

* `mcbo(algorithm = "BB")`: the branch-and-bound algorithm for the median
  constrained bucket order with `nbuckets` buckets did not always return the
  optimal solution. Three corrections:
  - a partial solution with k objects in nbk buckets is now discarded only if it
    cannot be completed with exactly `nbuckets` buckets, i.e. if
    nbk > nbuckets or nbk + (N - k) < nbuckets. The previous rule also
    discarded the partial solutions with fewer than `nbuckets` - 2 buckets,
    which could remove the optimal solution when it has a large bucket among
    the first objects inserted (e.g. N = 8, `nbuckets` = 4, optimal bucket
    sizes (5, 1, 1, 1): optimum missed in 35 out of 40 test data sets);
  - the threshold of the search is now the penalty of a feasible solution;
  - if no branch survives, the best known solution is returned.
  The corrected algorithm was checked against the exhaustive enumeration of all
  the bucket orders with `nbuckets` buckets (167 data sets with N = 6 to 9,
  with and without ties): optimal in all of them.

# ConsRank 3.1

## New features

* New function `bolt()` and new option `algorithm = "BOLT"` in `consrank()`:
  BOLT (Bucket-Order Local search with optimal Tie partitioning), a deterministic
  and scalable heuristic for the median ranking (Romano, D'Ambrosio and
  Conversano, manuscript under review). It works entirely on the combined input
  matrix: starting solutions computed from the matrix, a tie-aware insertion local
  search with moves evaluated in O(M) time, and the exact optimal tie partitioning
  of a fixed order in O(M^2) time (option `otp`; `otp = FALSE` gives the "fast"
  configuration of the paper).

  Details not specified in the paper, chosen in this implementation:
  - Markov chain score: random walk moving from each object to the objects that
    beat it, with weights the positive net preferences `max(c_ji - c_ij, 0)`; an
    object beaten by no other object gets a self-loop; a damping factor 0.85
    (uniform teleportation, as in PageRank) makes the stationary distribution
    unique. The damping factor is not in the paper.
  - Bucket pivot: the pivot is the object with the largest sum of absolute net
    preferences (least indecision, as in the paper); another object is tied with
    the pivot if tying is better than ordering them (`min(c_ij, c_ji) > 0`), is
    placed before it if it is preferred to it, and after it otherwise.
  - "Quick branch-and-bound start": the starting solution of the
    branch-and-bound algorithm (`findconsensusBB`), not the whole Quick
    algorithm, which is not feasible with a very large number of objects.
  - Net flow and Copeland starting solutions: objects with equal scores are tied.
  - No random perturbations of the starting solutions: in our experiments (12
    real and 180 simulated data sets) four perturbed starting solutions improved
    TauX in 8 cases out of 192, by at most 5.5e-5, at a cost of 10-15% more
    time, so neither `bolt()` nor `mcbo(algorithm = "BOLT")` uses them.
  - Local search: the objects are processed in the order of the current ranking;
    among equally good moves the first one is chosen; a move is accepted when it
    improves the objective by more than a small numerical tolerance.
  - Optimal tie partitioning of a ranking with ties: the objects of the same
    bucket are ordered by their index.

  Additions with respect to the paper:
  - `maxjudges`: maximum number of distinct rankings of the judges used as
    starting solutions, the most frequent first (default: all of them, as in the
    paper). Missing objects in these rankings are tied in the last position.
  - `full = TRUE`: search in the space of full rankings, with the local search
    without ties only (no tie partitioning).

  Because of these choices, when several optimal solutions exist the solutions
  returned may differ from those of the authors' code (e.g. on `EMD` two of the
  three optimal solutions are found, but not the same two reported in the paper).

* New arguments `cij` and `nj` in `consrank()`, `mcbo()` and `bolt()`: the
  combined input matrix (with the number of judges, or the sum of the weights,
  `nj`) can be given instead of the data `X`, since all the algorithms depend on
  the data only through this matrix. `X` now has default `NULL`; exactly one of
  `X` and `cij` must be given. With `cij`, `wk` is ignored (with a warning) and
  BOLT does not use the rankings of the judges as starting solutions. The results
  are identical to those obtained from the corresponding data.

* New option `algorithm = "BOLT"` in `mcbo()` (median constrained bucket order),
  not in D'Ambrosio et al. (2019): BOLT adapted to the constraint of exactly
  `nbuckets` buckets. For a fixed order, the optimal partition into exactly
  `nbuckets` buckets of consecutive objects is computed by dynamic programming in
  O(nbuckets M^2) time; from each starting solution (the same as in `bolt()`) the
  order is improved by the insertion local search without ties, then this optimal
  partition and a local search that keeps the number of buckets are alternated
  until no improvement. The local search has two families of moves: insertion
  moves that keep the number of buckets, then swaps of two objects of different
  buckets (the swap keeps every bucket size, so it is legal under any constraint
  on the buckets and finds improvements that need two objects to move together);
  the gain of a swap is computed in O(1) from cumulative sums over the buckets.
  Deterministic; option `maxjudges` as in `bolt()`. `bolt()` and
  `consrank(algorithm = "BOLT")` do not use the swap moves: their results are
  those of the BOLT algorithm of Romano, D'Ambrosio and Conversano.

* New argument `bucketsizes` in `mcbo()`: the consensus bucket order can be
  constrained by a pattern of bucket sizes instead of their number only. One
  element per bucket, a whole number greater than zero prescribing the size of
  that bucket and NA leaving it free; when the sizes given do not use all the
  objects and no NA is present, one free bucket is added at the end. The pattern
  contains as special cases the full rankings (a pattern of ones), the bucket
  orders with a given number of buckets (a pattern of NAs, which reproduces
  exactly `nbuckets` with `algorithm = "BOLT"`), the solution with all the
  objects tied (one bucket) and the top-k lists (k ones and a free bucket).
  It is handled by `algorithm = "BOLT"`, which is used by default when
  bucketsizes is given: for a fixed order of the objects the partition following
  the pattern is computed exactly by dynamic programming (one admissible cut
  where a size is prescribed, all of them where it is free), and the local search
  moves one object to another bucket making one object cross each boundary in
  between, so that every bucket keeps its size; when no such move improves, two
  objects of different buckets exchange their buckets (swap move), and the two
  families of moves are repeated until no swap improves. Either `nbuckets` or `bucketsizes`
  must be given, not both: to fix the number of buckets of a pattern, NA is used
  for the buckets of free size (e.g. `bucketsizes = c(1, 1, NA, NA, NA)`).

* New function `consranktopk()`: the median top-k list, i.e. the consensus
  ranking with k objects in strict order in the first k positions and all the
  other objects tied in the last position. It is a shortcut for
  `mcbo(..., bucketsizes = rep(1, k), algorithm = "BOLT")`, with the same
  results; the output also contains `TopK`, the names of the k objects in order.
  By default (`algorithm = "exact"`) the solution is exact: the agreement of a
  top-k list with the judges is a constant plus a linear term (column totals of
  the combined input matrix) plus the Kemeny agreement among the k objects of the
  list, so the objects that cannot belong to an optimal list are removed by two
  rules (exchange of an object of the list with one outside it; upper bound
  against a known list) and the best list of the remaining candidates is found by
  dynamic programming over their subsets with at most k objects; all the optimal
  lists are returned (at most 1000). When the candidates are more than `maxcand`
  (default 22; memory 8 bytes times 2^maxcand), the BOLT solution is returned.
  The output contains `Exact` (TRUE if certified optimal) and `Candidates`.
  `algorithm = "BOLT"` gives the BOLT solution only.

* New function `readpreflib()`: reads a data file of PrefLib (www.preflib.org),
  from a local file or a URL, and returns the rankings in the format of the
  package (one row per distinct ranking, one column per object named after the
  alternative, positions 1, 2, 3, ... without gaps) with their frequencies `wk`.
  It reads strict and tied orders, complete and incomplete (soc, soi, toc, toi),
  and categorical preferences (cat: the categories become the buckets, empty
  categories are skipped). The objects not ranked by a judge are tied in the last
  position (`unranked = "last"`, as in the files that PrefLib obtains by adding the
  unranked alternatives at the bottom) or left missing (`unranked = "na"`);
  rankings that become equal are merged. `expand = TRUE` gives one row per judge.

## Changes

* The notes "Combined Input Matrix contains only positive values: the median
  ranking is the all-tie solution" and "Combined Input Matrix contains only
  zeros: ..." of `consrank()` (all the algorithms) and `iwquickcons()` are now
  messages (`message()`) instead of printed values (`print()`): they can be
  suppressed with `suppressMessages()`, and they no longer appear as `[1] "..."`
  in the output of the functions that call `consrank()` (e.g. the trees of
  ConsRankClass). The results do not change.

* `mcbo()` accepts `nbuckets` equal to the number of objects: the solution is
  then a full ranking, and the result is the one of `consrank(..., full = TRUE)`
  with the same algorithm and parameters (before: error "nbuckets must be less
  than N").

* The defunct function `labels()` has been removed: being exported, it masked
  `base::labels()`, so that after `library(ConsRank)` any call to `labels()`
  stopped with the message that the function is defunct. Use `rank2order()`.

* `mcbo(algorithm = "quick")` is now deterministic, more accurate and faster;
  the results differ from those of the version described in D'Ambrosio et al.
  (2019). Changes: (1) the starting solution is the order given by the starting
  point of the branch-and-bound algorithm, cut into exactly `nbuckets` buckets
  by the optimal partition used by `algorithm = "BOLT"` (before: buckets of equal
  size by strength, plus three objects moved to random buckets from the second
  iteration on); (2) a candidate solution with a number of buckets different from
  `nbuckets` is brought back to `nbuckets` buckets by the optimal partition of its
  order (before: it was discarded, so that the search could hardly move); (3) the
  insertion passes are repeated until no improvement, `itermax` being the maximum
  number of passes (before: `itermax` independent single passes from perturbed
  starting solutions); (4) each solution is refined by the optimal partition of
  its order. The function `findconsensusBB_buckets()` (internal) has been removed.

* `mcbo(algorithm = "decor")` (differences with respect to D'Ambrosio et al., 2019):
  a good candidate is inserted in the initial population (the order of the
  starting point of the branch-and-bound algorithm, cut into exactly `nbuckets`
  buckets by the optimal partition), as in the DECoR algorithm of `consrank()`; a
  child with a wrong number of buckets is brought back to `nbuckets` buckets by the
  optimal partition of its order, instead of being replaced by a random
  individual; the final solution is refined by the optimal partition of its order.
  The differential evolution steps (mutation, crossover, selection) are unchanged.

* `kemenyd()` and `iw_kemenyd()` now stop with an error if the rankings contain
  missing values (NA): the Kemeny distance is not defined for partial rankings,
  which must be completed by the user (e.g. by tying the missing objects in the
  last position). Version 3.0 treated a missing object as tied with all the
  others, version 2.1.4 excluded the pairs with a missing object.

* `tau_x()` and `iw_tau_x()` with partial rankings (NA) now compute TauX as
  defined by Emond and Mason (2002): a pair of objects with a missing object has
  score 0 and does not contribute. Before, `tau_x()` was computed from the Kemeny
  distance, which codes a missing pair like a tie (so that, e.g., a pair missing
  in both rankings counted as an agreement), and `iw_tau_x()` ranked the missing
  objects last with distinct ranks. With complete rankings the results are
  unchanged.

* `scorematrix()` keeps the names of the objects also with the C++
  implementation (as the R implementation).

* Removed `src/Makevars.win`, which set the non-portable compilation flag `-O3`,
  and the unused duplicate C++ function `kemenydesign_cpp` (the R function
  `kemenydesign_cpp()` is unchanged).

* `polyplotnew()` has been removed: the improved version is now `polyplot()`.
  When `X` is given, the number of objects is the number of columns of `X` and
  `nobj` is ignored (with a warning if it is different); with more than 4 objects
  the function stops, because the polytope cannot be visualized.

* The package does not depend on `tidyr` any more. The only use of it was
  `tidyr::crossing()`, which built the table of all pairs of rankings in
  `iw_kemenyd()` and `iw_tau_x()` when two sets of rankings are compared; it is
  now built with base R. The results are unchanged. This removes `tidyr` and the
  packages it loads (including `vctrs`, `cli` and `rlang`) from the dependencies
  of ConsRank.

* The name of the algorithm of `consrank()` and `mcbo()` is now matched ignoring
  the case, so that `algorithm = "bolt"` and `algorithm = "BOLT"` are the same, and
  an unknown name gives a clear error message listing the available algorithms.
  Before, a name written with a different case, or a misspelled one, stopped with
  the message "object 'out' not found".

* `bolt()`, `consrank(algorithm = "BOLT")` and `mcbo(algorithm = "BOLT")` are faster
  and use less memory, with exactly the same results. Two changes: (1) the insertion
  local search rebuilt the order of the objects for every object of a sweep, and now
  rebuilds it only after an accepted move, because the solution does not change
  otherwise; (2) the matrices of the differences and of the sums of the combined input
  matrix are not built any more, so two n by n matrices less are allocated. The local
  search works directly on the combined input matrix, through the identity
  N(S) = constant - 2 * sum of c_ji over the ordered pairs of objects for which i is
  in an earlier bucket than j, the pairs of tied objects contributing only to the
  constant. On simulated data with 20 judges, the time goes from 3.4 to 0.9 seconds
  with 500 objects and from 14.6 to 4.0 seconds with 1000 objects.

## Bug fixes

* `iwquickcons()` with a single judge (a matrix with one row) stopped with an
  error ("dim(X) must have a positive length" when `w` is a vector; "object
  'cij' not found" otherwise). Now the median ranking is the ranking of the
  judge (with `full = TRUE` and ties, the ties are broken in the order of the
  columns, as in `consrank()`), with its item-weighted TauX. Results with two
  or more judges do not change.
* `iwquickcons()` with `full = TRUE` could return a ranking with ties: the
  starting solution of the search was computed with ties (without the option
  `FULL`, which `consrank(algorithm = "quick")` uses) and the search kept them.
  Now the starting solution has no ties, as in `consrank()`, and the solution is
  a full ranking (on Univranks with item weights, the same TauX as the best of
  all the 720 permutations). With `full = FALSE` the results do not change.

* `consrank()` with a single judge: with `algorithm = "BB"` and `full = TRUE`
  the TauX returned was 0 instead of 1; with `algorithm = "decor"` the function
  stopped with the error "object 'Consensus' not found". With `full = TRUE` and
  a judge's ranking with ties, all the algorithms (also `bolt()`) returned that
  ranking, which is not a full ranking: now the ties are broken in the order of
  the columns (every full ranking compatible with the judge's ranking is a
  median ranking) and its TauX is returned. Results with two or more judges do
  not change.
* `mcbo(algorithm = "decor")`: when several individuals of the initial
  population had the same best cost, the stored best individuals were shifted
  with respect to their costs, so that the solution returned could differ from
  the one with the reported TauX.
* `consrank(algorithm = "decor")` stopped with an error when several individuals of
  the initial population had the same best cost.
* `consrank(algorithm = "decor")`: when the combined input matrix contains only
  zeros the consensus is now a 1 x M matrix with the names of the objects, as for
  the other algorithms (before: a vector).
* `consrank(algorithm = "fast")`: wrong condition `nrow(CR>1)` (now `nrow(CR) > 1`);
  the results do not change.

* Non-integer weights (e.g. `wk / sum(wk)` instead of `wk`): the branch-and-bound
  and Quick algorithms of `consrank()` and `mcbo()` compared the penalties exactly,
  and used the sign of elements of the combined input matrix that are zero up to
  rounding errors. They could therefore stop with an error ("argument of length
  0") or return a different solution than with the equivalent integer weights.
  The comparisons now use a tolerance relative to the scale of the combined input
  matrix, and its elements that are zero up to rounding errors are set to zero.
  With integer weights the results are unchanged; with weights multiplied by any
  constant the results are now the same as with the original weights.

* `mcbo(algorithm = "quick")` could return a solution with a number of buckets
  different from `nbuckets` (18 times out of 100 in a test with 4 to 7 objects).

* `polyplot()`: the point B > A > C of the polytope of 3 objects was missing
  (the ranking BAC was coded as (3, 1, 3)), and the weights `Wk` were assigned to
  the wrong points. The area of each symbol is now proportional to the frequency
  of the ranking (argument `size`), with semi-transparent symbols.

* `univranks()` stopped with an error when the input defined a single ranking or a
  single partition: a number of items equal to one (`univranks(1)`), a number of
  buckets equal to one or to the number of items (e.g. `univranks(4, 1)`,
  `univranks(4, 4)`), or a matrix with one row. `partitions(1, 1)` returned the
  ordering "{a}" instead of "a". The other results are unchanged. A numeric
  matrix was taken as a matrix of orderings, so that a ranking with ties given as a
  one-row matrix stopped with an error: a one-row numeric matrix is now the
  corresponding rank vector, and a numeric matrix with more rows gives a clear error
  (the input is a number of items, a single ranking or an ordering).

* Defunct functions: the messages of `EMCons()` and `FASTDECOR()` pointed to
  algorithms that do not exist (`"EM"`, `"fastdecor"`); they now point to
  `consrank(algorithm = "BB")` and `consrank(algorithm = "decor")`.

* Help page of the package (`?ConsRank`): removed the version and date, which
  were those of version 2.1.0, and completed the list of contributors as in the
  DESCRIPTION file.
