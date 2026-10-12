#' Median top-k list
#'
#' Find the median top-k list according to the Kemeny's axiomatic approach: the
#' consensus ranking in which k objects are ranked in strict order in the first k
#' positions and all the other objects are tied in the last position. It is the
#' median constrained bucket order with bucket sizes 1, ..., 1 (k times) and a
#' last bucket with the remaining objects (see \code{\link{mcbo}} with
#' \code{bucketsizes = rep(1, k)}). By default the solution is exact whenever
#' possible (see Details).
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to
#'   be judged. Each row is a ranking of the objects which are represented by
#'   the columns. If X contains the rankings observed only once, the argument wk
#'   can be used. Either X or cij must be given
#' @param k The number of objects of the top-k list: a whole number between 1 and
#'   M-1 (with k = M-1 the top-k list is a full ranking)
#' @param wk Optional: the frequency of each ranking in the data
#' @param ps Logical. If TRUE, displays information on screen. Default TRUE
#' @param maxjudges Maximum number of distinct rankings of the judges used as
#'   starting solutions by BOLT, the most frequent first. Default Inf (all of
#'   them); 0 uses only the starting solutions computed from the combined input
#'   matrix
#' @param cij The combined input matrix (M by M, see \code{\link{combinpmatr}}),
#'   to be given instead of X. With cij, the argument wk is ignored and the
#'   rankings of the judges are not used as starting solutions
#' @param nj The number of judges (the sum of the weights, if any) that produced
#'   cij. Required when cij is given: it is used to compute the averaged TauX
#' @param algorithm \code{"exact"} (default): the exact solution, found by
#'   dynamic programming on the objects that can belong to an optimal list, when
#'   they are at most maxcand; otherwise the BOLT solution. \code{"BOLT"}: always
#'   the BOLT solution, as \code{mcbo(X, bucketsizes = rep(1, k), algorithm =
#'   "BOLT")}
#' @param maxcand The largest number of candidate objects for the exact
#'   solution (default 22). The memory used grows as 2^maxcand (8 bytes times
#'   2^maxcand: 32 MB with 22); at most 30
#'
#' @return A list containing the following components:
#' \tabular{lll}{
#'   TopK \tab \tab the k objects of the median top-k list, in order (a matrix
#'   with one row per solution if there are several)\cr
#'   Consensus \tab \tab the consensus ranking(s) of all the objects, with the
#'   objects outside the top-k list tied in position k+1\cr
#'   Tau \tab \tab averaged TauX rank correlation coefficient\cr
#'   Exact \tab \tab TRUE if the solution is certified optimal, FALSE if it is
#'   the BOLT solution\cr
#'   Candidates \tab \tab the number of objects that can belong to an optimal
#'   list (NA with algorithm = "BOLT")\cr
#'   Eltime \tab \tab elapsed time in seconds
#' }
#'
#' @details
#' The objects of the top-k list are those for which the judges agree most on
#' placing them in front of the others, and their order is the one that best
#' agrees with the judges; the objects outside the list are tied, so that their
#' mutual order does not matter. Like the other functions of the package, the
#' solution maximizes the averaged TauX rank correlation coefficient with the
#' rankings of the judges. The procedure is deterministic.
#'
#' With the combined input matrix C (elements c_ij), the agreement of a top-k
#' list with the judges is a constant plus
#' \deqn{\sum_{t=1}^{k} g_{a_t} + 2 \sum_{s<t} c_{a_s a_t},}
#' where \eqn{a_1, \ldots, a_k} are the objects of the list in order and
#' \eqn{g_a = -2 \sum_{j \neq a} c_{ja}} depends only on the column total of
#' the object a. With \code{algorithm = "exact"}, the objects that cannot belong
#' to an optimal list are first removed, by two rules applied repeatedly: an
#' object of an optimal list cannot be improved by replacing it with an object
#' outside the list, and the value of a list containing the object, bounded from
#' above, must reach the value of a known list. The best list of the remaining
#' candidates is then found exactly by dynamic programming over their subsets with
#' at most k objects, and all the optimal lists are returned (at most 1000). The
#' first rule needs no known list; if the candidates are more than maxcand, the
#' BOLT solution is computed and used in the second rule; if they are still more
#' than maxcand, the BOLT solution is returned and Exact is FALSE.
#'
#' @examples
#' data(EMD)
#' CR <- consranktopk(EMD[, 1:15], k = 3, wk = EMD[, 16])
#' CR$TopK
#' CR$Tau
#' CR$Exact
#'
#' @references
#' Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with
#' application to the consensus ranking problem. Journal of Multi-Criteria
#' Decision Analysis, 11(1), 17-28.
#'
#' D'Ambrosio, A., Iorio, C., Staiano, M., and Siciliano, R. (2019).
#' Median constrained bucket order rank aggregation.
#' Computational Statistics, 34(2), 787-802.
#' \doi{10.1007/s00180-018-0858-z}
#'
#' Romano, M., D'Ambrosio, A., and Conversano, C. A deterministic and scalable
#' heuristic for median ranking detection: tie-aware local search with optimal
#' tie partitioning on the combined input matrix. Manuscript under review.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{mcbo}} for bucket orders with a given number of buckets or
#'   a given pattern of bucket sizes
#' @seealso \code{\link{consrank}} for the unconstrained median ranking
#' @seealso \code{\link{bolt}} for the BOLT algorithm
#'
#' @keywords Median ranking
#' @keywords Top-k list
#'
#' @export

consranktopk <- function(X = NULL, k, wk = NULL, ps = TRUE, maxjudges = Inf,
                         cij = NULL, nj = NULL, algorithm = c("exact", "BOLT"),
                         maxcand = 22) {

  # the algorithm: "exact" by default (case insensitive, as in consrank)
  algorithm <- matchalgorithm(algorithm[1], c("exact", "BOLT"))
  # number of objects, from the data or from the combined input matrix
  M <- if (!is.null(X)) NCOL(X) else NCOL(cij)
  # k must be one whole number
  if (missing(k) || length(k) != 1 || !is.numeric(k) || is.na(k) || k != round(k)) {
    stop("k must be a whole number: the number of objects of the top-k list")
  }
  # k must leave at least one object outside the list
  if (k < 1 || k > M - 1) {
    stop("k must be between 1 and the number of objects minus one (", M - 1, "); ",
         "for a full ranking use consrank(..., full = TRUE)")
  }
  # maxcand must be a whole number between 1 and 30 (memory: 8 bytes x 2^maxcand)
  if (length(maxcand) != 1 || !is.numeric(maxcand) || is.na(maxcand) ||
      maxcand != round(maxcand) || maxcand < 1 || maxcand > 30) {
    stop("maxcand must be a whole number between 1 and 30: the exact solution ",
         "needs 8 bytes times 2^maxcand of memory")
  }

  # the BOLT solution: the bucket order with k buckets of one object, the free
  # bucket with the remaining objects added at the end by mcbo
  runbolt <- function(quiet) {
    mcbo(X, wk = wk, ps = ps && !quiet, algorithm = "BOLT", maxjudges = maxjudges,
         cij = cij, nj = nj, bucketsizes = rep(1, k))
  }

  # the consensus as a matrix with one row per solution
  asmatrix <- function(cons) {
    if (is.null(dim(cons))) cons <- matrix(cons, nrow = 1, dimnames = list(NULL, names(cons)))
    cons
  }

  # the top-k lists of a matrix of consensus rankings (one row per solution)
  topklists <- function(cons) {
    # the names of the objects, or their column numbers if there are no names
    labs <- if (is.null(colnames(cons))) as.character(seq_len(M)) else colnames(cons)
    # for each solution, the objects in the first k positions, in order
    topk <- t(apply(cons, 1, function(r) labs[order(r)[seq_len(k)]]))
    # apply returns a vector of length k when k = 1: back to one column
    if (k == 1) topk <- matrix(topk, ncol = 1)
    # the columns are the positions 1, ..., k
    colnames(topk) <- paste0("pos", seq_len(k))
    # one solution: a plain vector of names
    if (nrow(topk) == 1) topk <- topk[1, ]
    topk
  }

  # the output of the BOLT solution: the top-k list first, then the output of
  # mcbo (Exact = FALSE, Candidates = NA)
  boltoutput <- function(out) {
    cons <- asmatrix(out$Consensus)
    c(list(TopK = topklists(cons)), out[c("Consensus", "Tau")],
      list(Exact = FALSE, Candidates = NA), out["Eltime"])
  }

  # ---- algorithm = "BOLT": the BOLT solution only ----
  if (algorithm == "BOLT") return(boltoutput(runbolt(quiet = FALSE)))

  # ---- algorithm = "exact" ----
  tic <- proc.time()[3]
  # the combined input matrix and the number of judges
  cijgiven <- !is.null(checkXcij(X, cij, nj, wk))
  if (cijgiven) {
    C <- as.matrix(cij)
    njudges <- nj
    labs <- colnames(C)
  } else {
    # the data as a matrix with one row per ranking
    if (is(X, "data.frame")) X <- as.matrix(X)
    if (is.null(nrow(X))) X <- matrix(X, nrow = 1)
    # as in mcbo, at least two rankings
    if (nrow(X) < 2) stop("Data matrix must contain at least two different rankings")
    # the weights must be one positive number per ranking
    if (!is.null(wk)) {
      if (length(wk) != nrow(X)) stop("Length of wk must equal number of rows in X")
      if (any(wk <= 0)) stop("All weights must be positive")
    }
    C <- combinpmatr(X, Wk = wk)
    njudges <- if (is.null(wk)) nrow(X) else sum(wk)
    labs <- colnames(X)
  }
  # elements equal to zero up to rounding errors
  C <- cleancij(C)
  # no information at all (all elements zero): the BOLT solution, as in mcbo
  if (all(C == 0)) return(boltoutput(runbolt(quiet = FALSE)))

  # numerical tolerance, as in BOLT
  tol <- 1e-10 * max(1, sum(abs(C)))
  # the agreement of the ranking with all the objects tied (constant term)
  tiedval <- sum((C + t(C))[upper.tri(C)])

  # first rule only: the objects that can belong to an optimal list
  cand <- topk_cand_impl(C, as.integer(k), -Inf, tol)$cand
  bolt <- NULL
  # too many candidates: the value of the BOLT solution gives the second rule
  if (sum(cand) > maxcand) {
    bolt <- runbolt(quiet = TRUE)
    bcons <- asmatrix(bolt$Consensus)
    # value of the BOLT list (agreement minus the constant term)
    F0 <- bolt_obj_impl(t(bcons[1, , drop = FALSE]), C) - tiedval
    cand <- topk_cand_impl(C, as.integer(k), F0, tol)$cand
  }
  # number of candidates
  ncand <- sum(cand)

  if (ncand <= maxcand) {
    # exact solution: dynamic programming on the candidates (numbered from 0)
    dp <- topk_dp_impl(C, as.integer(which(cand) - 1L), as.integer(k), tol, 1000L)
    if (dp$truncated) {
      warning("there are more than 1000 optimal top-k lists: only 1000 are returned")
    }
    # one consensus ranking per optimal list: the list in positions 1, ..., k,
    # the other objects tied in position k + 1
    seqs <- dp$sequences
    cons <- matrix(k + 1, nrow = nrow(seqs), ncol = M, dimnames = list(NULL, labs))
    for (i in seq_len(nrow(seqs))) cons[i, seqs[i, ]] <- seq_len(k)
    exact <- TRUE
  } else {
    # still too many candidates: the BOLT solution, not certified
    cons <- asmatrix(bolt$Consensus)
    exact <- FALSE
  }

  # averaged TauX of each solution (they are all equal)
  tau <- bolt_obj_impl(t(cons), C) / (njudges * M * (M - 1))
  toc <- proc.time()[3]

  # information on screen
  if (ps) {
    if (exact) {
      message("exact top-k: ", ncand, " candidate objects, ", nrow(cons),
              " optimal solution(s) with Tau = ", round(tau[1], 5))
    } else {
      message("exact top-k not possible: ", ncand, " candidate objects (more than ",
              "maxcand = ", maxcand, "); BOLT solution with Tau = ", round(tau[1], 5))
    }
  }

  list(TopK = topklists(cons), Consensus = cons,
       Tau = if (nrow(cons) == 1) tau else matrix(tau, ncol = 1),
       Exact = exact, Candidates = ncand, Eltime = toc - tic)
}
