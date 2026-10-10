#' Median top-k list
#'
#' Find the median top-k list according to the Kemeny's axiomatic approach: the
#' consensus ranking in which k objects are ranked in strict order in the first k
#' positions and all the other objects are tied in the last position. It is the
#' median constrained bucket order with bucket sizes 1, ..., 1 (k times) and a
#' last bucket with the remaining objects, found by the BOLT algorithm (see
#' \code{\link{mcbo}} with \code{bucketsizes = rep(1, k)}).
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
#'   starting solutions, the most frequent first. Default Inf (all of them);
#'   0 uses only the starting solutions computed from the combined input matrix
#' @param cij The combined input matrix (M by M, see \code{\link{combinpmatr}}),
#'   to be given instead of X. With cij, the argument wk is ignored and the
#'   rankings of the judges are not used as starting solutions
#' @param nj The number of judges (the sum of the weights, if any) that produced
#'   cij. Required when cij is given: it is used to compute the averaged TauX
#'
#' @return A list containing the following components:
#' \tabular{lll}{
#'   TopK \tab \tab the k objects of the median top-k list, in order (a matrix
#'   with one row per solution if there are several)\cr
#'   Consensus \tab \tab the consensus ranking(s) of all the objects, with the
#'   objects outside the top-k list tied in position k+1\cr
#'   Tau \tab \tab averaged TauX rank correlation coefficient\cr
#'   Eltime \tab \tab elapsed time in seconds
#' }
#'
#' @details
#' The function is a shortcut for \code{mcbo(X, wk = wk, bucketsizes = rep(1, k),
#' algorithm = "BOLT")}: the results are the same. The objects of the top-k list
#' are those for which the judges agree most on placing them in front of the
#' others, and their order is the one that best agrees with the judges; the
#' objects outside the list are tied, so that their mutual order does not matter.
#' Like the other functions of the package, the solution maximizes the averaged
#' TauX rank correlation coefficient with the rankings of the judges. The
#' procedure is deterministic.
#'
#' @examples
#' data(EMD)
#' CR <- consranktopk(EMD[, 1:15], k = 3, wk = EMD[, 16])
#' CR$TopK
#' CR$Tau
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
                         cij = NULL, nj = NULL) {

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

  # the top-k list is the bucket order with k buckets of one object; the free
  # bucket with the remaining objects is added at the end by mcbo
  out <- mcbo(X, wk = wk, ps = ps, algorithm = "BOLT", maxjudges = maxjudges,
              cij = cij, nj = nj, bucketsizes = rep(1, k))

  # the consensus as a matrix with one row per solution
  cons <- out$Consensus
  if (is.null(dim(cons))) cons <- matrix(cons, nrow = 1, dimnames = list(NULL, names(cons)))
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

  # the top-k list first, then the output of mcbo
  c(list(TopK = topk), out)
}
