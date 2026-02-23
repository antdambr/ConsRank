#' Median Constrained Bucket Order (MCBO)
#'
#' Find the median ranking constrained to exactly b buckets (tied groups) according 
#' to Kemeny's axiomatic approach. This implements the algorithms described in 
#' D'Ambrosio et al. (2019) for rank aggregation with fixed bucket constraints.
#'
#' @param X A N by M data matrix, in which there are N judges and M objects to be judged. 
#'   Each row is a ranking of the objects which are represented by the columns. 
#'   If X contains the rankings observed only once, the argument wk can be used
#' @param nbuckets Integer. The number of buckets (tied groups) the consensus ranking 
#'   must contain. Must be between 2 and M-1 (where M is the number of objects)
#' @param wk Optional: the frequency of each ranking in the data
#' @param ps Logical. If TRUE, displays progress information on screen. Default TRUE
#' @param algorithm Character string specifying the algorithm to use. One of:
#'   \itemize{
#'     \item "BB" - Branch-and-Bound (exact, default). Best for M <= 15
#'     \item "quick" - Quick algorithm (heuristic). Best for 15 < M <= 50
#'     \item "decor" - Differential Evolution (metaheuristic). Best for M > 50
#'   }
#' @param itermax Integer. Maximum number of iterations for "quick" and "decor" algorithms. 
#'   Default 10
#' @param np Integer. Number of population individuals for "decor" algorithm. Default 10
#' @param gl Integer. Generations limit for "decor": maximum number of consecutive 
#'   generations without improvement. Default 100
#' @param ff Numeric. Scaling rate for mutation in "decor". Must be in [0,1]. Default 0.4
#' @param cr Numeric. Crossover range for "decor". Must be in [0,1]. Default 0.8
#' @param use_cpp Logical. If TRUE (default), use optimized C++ implementations 
#'   for core functions (combinpmatr, scorematrix, PenaltyBB2, etc.)
#'
#' @return A list containing the following components:
#' \tabular{lll}{
#'   Consensus \tab \tab The consensus ranking with exactly nbuckets buckets\cr
#'   Tau \tab \tab Averaged TauX rank correlation coefficient\cr
#'   Eltime \tab \tab Elapsed time in seconds
#' }
#'
#' @details 
#' The median constrained bucket order problem finds the ranking that minimizes the 
#' sum of Kemeny distances to the input rankings, subject to the constraint that 
#' the solution must have exactly \code{nbuckets} groups of tied items.
#' 
#' This is useful in applications where the output must conform to predetermined 
#' categories, such as:
#' \itemize{
#'   \item Wine quality classifications (e.g., 5 fixed tiers)
#'   \item Medical triage systems (fixed severity codes)
#'   \item Educational grading (fixed letter grades: A, B, C, D, F)
#' }
#' 
#' The search space is restricted to Z^{n\textbackslash b}, which contains
#' \deqn{\sum_{b=1}^{n} b! \cdot S(n,b)}
#' rankings, where S(n,b) is the Stirling number of the second kind.
#' 
#' \strong{Algorithm Selection Guidelines:}
#' \itemize{
#'   \item \strong{BB}: Exact solution, guaranteed optimal. Use for M <= 15 items
#'   \item \strong{quick}: Fast heuristic, near-optimal. Use for 15 < M <= 50 items
#'   \item \strong{decor}: Best for large problems. Use for M > 50 items
#' }
#' 
#' For stochastic algorithms (quick, decor), consider running multiple times 
#' (controlled by \code{itermax}) to avoid local optima.
#'
#' @examples
#' \dontrun{
#' # Simple example with 10 judges ranking 5 items into 3 buckets
#' data(Idea)
#' RevIdea <- 6 - Idea  # Reverse ranking
#' CR <- mcbo(RevIdea, nbuckets = 3, algorithm = "BB")
#' print(CR$Consensus)
#' print(CR$Tau)
#' 
#' # Larger dataset with Quick algorithm
#' data(EMD)
#' CR_quick <- mcbo(EMD[,1:15], nbuckets = 5, wk = EMD[,16], 
#'                  algorithm = "quick", itermax = 20)
#' 
#' # Very large dataset with DECoR
#' # CR_decor <- mcbo(large_data, nbuckets = 7, algorithm = "decor",
#' #                  np = 20, itermax = 100)
#' }
#'
#' @references 
#' D'Ambrosio, A., Iorio, C., Staiano, M., and Siciliano, R. (2019). 
#' Median constrained bucket order rank aggregation. 
#' Computational Statistics, 34(2), 787-802. 
#' \doi{10.1007/s00180-018-0858-z}
#' 
#' Emond, E. J., and Mason, D. W. (2002). 
#' A new rank correlation coefficient with application to the consensus ranking problem. 
#' Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' Kemeny, J. G., and Snell, L. J. (1962). 
#' Preference ranking: an axiomatic approach. 
#' Mathematical models in the social sciences, 9-23.
#'
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#'
#' @seealso \code{\link{consrank}} for unconstrained consensus ranking
#' @seealso \code{\link{combinpmatr}} for computing the combined input matrix
#' @seealso \code{\link{scorematrix}} for computing score matrices
#' @seealso \code{\link{tau_x}} for TauX correlation coefficient
#' @seealso \code{\link{kemenyd}} for Kemeny distance
#' @seealso \code{\link{stirling2}} for Stirling numbers (bucket combinatorics)
#'
#' @keywords Median ranking
#' @keywords Consensus ranking
#' @keywords Bucket order
#' @keywords Kemeny distance
#' @keywords Rank aggregation
#'
#' @export

mcbo <- function(X, nbuckets, wk = NULL, ps = TRUE, 
                 algorithm = "BB", itermax = 10,
                 np = 10, gl = 100, ff = 0.4, cr = 0.8,
                 use_cpp = TRUE) {
  
  # ══════════════════════════════════════════════════════════════════════════
  # 1. INPUT VALIDATION AND STANDARDIZATION
  # ══════════════════════════════════════════════════════════════════════════
  
  # Convert data.frame to matrix
  if (is(X, "data.frame")) {
    X <- as.matrix(X)
  }
  
  # Handle vector input (single ranking)
  if (is(nrow(X), "NULL")) {
    X <- matrix(X, nrow = 1)
  }
  
  M <- nrow(X)  # Number of judges
  N <- ncol(X)  # Number of objects
  
  # ══════════════════════════════════════════════════════════════════════════
  # 2. VALIDATION CHECKS
  # ══════════════════════════════════════════════════════════════════════════
  
  # Check: Need at least 2 rankings
  if (M < 2) {
    stop("Data matrix must contain at least two different rankings")
  }
  
  # Check: nbuckets must be integer
  if (!is.numeric(nbuckets) || length(nbuckets) != 1 || nbuckets != as.integer(nbuckets)) {
    stop("nbuckets must be a single integer value")
  }
  nbuckets <- as.integer(nbuckets)
  
  # Check: nbuckets range
  if (nbuckets < 2) {
    stop("nbuckets must be at least 2 (for nbuckets=1, the trivial solution is all items tied)")
  }
  
  if (nbuckets >= N) {
    stop(paste0("nbuckets must be less than N (number of objects = ", N, ").\n",
                "For nbuckets = N, use consrank() with full = TRUE"))
  }
  
  # Check: algorithm validity
  valid_algorithms <- c("BB", "quick", "decor")
  if (!algorithm %in% valid_algorithms) {
    stop(paste0("algorithm must be one of: ", paste(valid_algorithms, collapse = ", ")))
  }
  
  # Check: weights validity
  if (!is(wk, "NULL")) {
    if (is.numeric(wk)) {
      wk <- matrix(wk, ncol = 1)
    }
    if (length(wk) != M) {
      stop("Length of wk must equal the number of rows in X")
    }
    if (any(wk <= 0)) {
      stop("All weights in wk must be positive")
    }
  }
  
  # Check: parameter ranges for DECOR
  if (algorithm == "decor") {
    if (ff < 0 || ff > 1) {
      stop("ff (mutation scaling rate) must be in [0, 1]")
    }
    if (cr < 0 || cr > 1) {
      stop("cr (crossover range) must be in [0, 1]")
    }
    if (np < 1) {
      stop("np (population size) must be at least 1")
    }
    if (gl < 1) {
      stop("gl (generations limit) must be at least 1")
    }
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # 3. COMPUTE COMBINED INPUT MATRIX
  # ══════════════════════════════════════════════════════════════════════════
  
  cij <- combinpmatr(X, Wk = wk, use_cpp = use_cpp)
  
  # Compute number of judges (for TauX calculation)
  nj <- if (is(wk, "NULL")) M else sum(wk)
  
  # ══════════════════════════════════════════════════════════════════════════
  # 4. CHECK FOR TRIVIAL SOLUTIONS
  # ══════════════════════════════════════════════════════════════════════════
  
  # Check 1: All zeros in CIJ → any ranking is optimal
  if (sum(cij == 0) == length(cij)) {
    message("Combined Input Matrix contains only zeros: any ranking with ", 
            nbuckets, " buckets is a median ranking")
    
    # Generate a simple solution: assign items round-robin to buckets
    cr <- matrix(rep(1:nbuckets, length.out = N), nrow = 1)
    colnames(cr) <- colnames(X)
    
    return(list(Consensus = cr, Tau = 0, Eltime = 0))
  }
  
  # Check 2: All positive in CIJ (for non-bucket case, would give all-tie)
  # For bucket case, we still need to compute optimal assignment
  # So we don't return early here
  
  # ══════════════════════════════════════════════════════════════════════════
  # 5. ALGORITHM DISPATCH
  # ══════════════════════════════════════════════════════════════════════════
  
  tic <- proc.time()[3]
  
  if (algorithm == "BB") {
    
    # Branch-and-Bound: exact solution
    out <- BB_buckets(cij = cij, 
                      nbuckets = nbuckets, 
                      nj = nj, 
                      Wk = wk, 
                      PS = ps,
                      use_cpp = use_cpp)
    
  } else if (algorithm == "quick") {
    
    # Quick algorithm: fast heuristic
    out <- Quick_buckets(cij = cij, 
                         nbuckets = nbuckets, 
                         nj = nj, 
                         Wk = wk, 
                         maxiter = itermax,
                         PS = ps,
                         use_cpp = use_cpp)
    
  } else if (algorithm == "decor") {
    
    # Differential Evolution: metaheuristic for large problems
    out <- DECoR_buckets(cij = cij, 
                         nbuckets = nbuckets, 
                         nj = nj, 
                         Wk = wk, 
                         maxiter = itermax,
                         NP = np, 
                         L = gl, 
                         FF = ff, 
                         CR = cr,
                         PS = ps,
                         use_cpp = use_cpp)
    
  }
  
  toc <- proc.time()[3]
  
  # ══════════════════════════════════════════════════════════════════════════
  # 6. FINALIZE OUTPUT
  # ══════════════════════════════════════════════════════════════════════════
  
  out$Eltime <- toc - tic
  
  # Ensure column names are preserved
  if (!is.null(colnames(X))) {
    colnames(out$Consensus) <- colnames(X)
  }
  
  return(out)
}

#-------------------------------------------
#Branch-and-Bound algorithm for Median Constrained Bucket Order
#
# Exact algorithm to find the median ranking with exactly nbuckets buckets.
# This is the core implementation for the BB algorithm in mcbo().


BB_buckets <- function(cij, nbuckets, nj, Wk = NULL, PS = TRUE, use_cpp = TRUE) {
  
  N <- ncol(cij)
  isw <- is.null(Wk)  # TRUE if no weights
  
  tic <- proc.time()[3]
  
  # ══════════════════════════════════════════════════════════════════════════
  # STEP 1: Generate initial candidate using quantile-based heuristic
  # ══════════════════════════════════════════════════════════════════════════
  
  if (PS) {
    message("Generating initial candidate...")
  }
  
  # Generate single deterministic initial candidate
  R_init <- findconsensusBB_buckets(cij, nbuckets)
  
  Sij_init <- scorematrix(R_init, use_cpp = use_cpp)
  Po_init <- sum(abs(cij)) - sum(cij * Sij_init)  
  
  if (PS) {
    message("Initial candidate buckets: ", length(unique(R_init[1,])))
    message("Starting Branch-and-Bound refinement...")
  }
  
  # ══════════════════════════════════════════════════════════════════════════
  # STEP 2: Refine with full Branch-and-Bound
  # ══════════════════════════════════════════════════════════════════════════
  
  consensus <- BBconsensus_buckets(
    RR = R_init, 
    cij = cij, 
    Po = Po_init, 
    nbuckets = nbuckets,
    PS = PS,
    use_cpp = use_cpp
  )
  
  # ══════════════════════════════════════════════════════════════════════════
  # STEP 3: Compute TauX correlation coefficient
  # ══════════════════════════════════════════════════════════════════════════
  
  if (nrow(consensus) == 1) {
    
    # Single solution
    Sij <- scorematrix(consensus, use_cpp = use_cpp)
    
    if (!isw) {
      TauX <- sum(cij * Sij) / (sum(Wk) * (N * (N - 1)))
    } else {
      TauX <- sum(cij * Sij) / (nj * (N * (N - 1)))
    }
    
  } else {
    
    # Multiple solutions
    TauX <- numeric(nrow(consensus))
    
    for (k in 1:nrow(consensus)) {
      Sij <- scorematrix(matrix(consensus[k, ], nrow = 1), use_cpp = use_cpp)
      
      if (!isw) {
        TauX[k] <- sum(cij * Sij) / (sum(Wk) * (N * (N - 1)))
      } else {
        TauX[k] <- sum(cij * Sij) / (nj * (N * (N - 1)))
      }
    }
  }
  
  toc <- proc.time()[3]
  eltime <- toc - tic
  
  return(list(
    Consensus = reordering(consensus), 
    Tau = TauX, 
    Eltime = eltime
  ))
}


# ══════════════════════════════════════════════════════════════════════════
# HELPER FUNCTION: Quick initialization for BB
# ══════════════════════════════════════════════════════════════════════════

findconsensusBB_buckets <- function(cij, nbuckets) {
  
  N <- ncol(cij)
  
  # Forza di ogni item
  strength <- colSums(cij)
  
  # Dividi in nbuckets usando quantili
  breaks <- quantile(strength, probs = seq(0, 1, length.out = nbuckets + 1))
  
  # Assegna bucket
  bucket <- cut(strength, breaks = breaks, labels = FALSE, 
                include.lowest = TRUE)
  
  # Fix se mancano buckets
  actual_buckets <- length(unique(bucket))
  
  if (actual_buckets < nbuckets) {
    # Splitta i bucket più grandi fino ad avere nbuckets
    while (length(unique(bucket)) < nbuckets) {
      bucket_sizes <- table(bucket)
      largest_bucket <- as.numeric(names(which.max(bucket_sizes)))
      items_in_largest <- which(bucket == largest_bucket)
      
      # Splitta a metà
      mid <- length(items_in_largest) %/% 2
      new_bucket_id <- max(bucket) + 1
      bucket[items_in_largest[(mid+1):length(items_in_largest)]] <- new_bucket_id
    }
  } else if (actual_buckets > nbuckets) {
    # Merge i bucket più piccoli
    while (length(unique(bucket)) > nbuckets) {
      bucket_sizes <- table(bucket)
      smallest_bucket <- as.numeric(names(which.min(bucket_sizes)))
      items_in_smallest <- which(bucket == smallest_bucket)
      
      # Trova bucket adiacente per merge
      bucket_means <- tapply(strength, bucket, mean)
      adjacent <- which.min(abs(bucket_means - bucket_means[smallest_bucket]))
      
      bucket[items_in_smallest] <- adjacent
    }
  }
  
  # Riordina bucket per strength media
  bucket_means <- tapply(strength, bucket, mean)
  bucket_order <- rank(bucket_means)
  
  final_ranking <- bucket_order[bucket]
  
  return(matrix(final_ranking, 1, N))
}


# ══════════════════════════════════════════════════════════════════════════
# BBconsensus_buckets: Find BEST solution with EXACTLY nbuckets
# (NOT the globally optimal solution - that's what consrank() does)
# ══════════════════════════════════════════════════════════════════════════

BBconsensus_buckets <- function(RR, cij, Po, nbuckets, PS = TRUE, use_cpp = TRUE) {
  
  CR <- RR
  a <- t(matrix(sort(RR, decreasing = TRUE)))
  ord <- t(matrix(order(RR, decreasing = TRUE)))
  r <- ReorderingBB(RR, use_cpp = use_cpp)
  
  BR.R <- r       
  BR.P <- 0       
  WCi <- 1        
  lambda <- 1     
  nobj <- ncol(RR)
  
  # ═══════════════════════════════════════════════════════════════════════
  # CONVERGENCE: Track best k-bucket solution found
  # ═══════════════════════════════════════════════════════════════════════
  best_kbucket_penalty <- Inf
  best_kbucket_solution <- NULL
  no_improvement_rounds <- 0
  
  while (WCi == 1) {
    
    if (PS) {
      message("BB Round ", lambda)
    }
    
    # ═══════════════════════════════════════════════════════════════════════
    # PRIMARY LOOP: Add k-th best ranked object
    # ═══════════════════════════════════════════════════════════════════════
    
    for (k in 2:ncol(a)) {
      
      B <- nrow(BR.R)
      b <- 1:k
      
      KR.R_list <- vector("list", B)
      KR.P_list <- vector("list", B)
      valid_count <- 0
      
      for (nb in 1:B) {
        BR.R[nb, ] <- ReorderingBB(matrix(BR.R[nb, ], nrow = 1), use_cpp = use_cpp)
        
        rpbr <- branches_buckets(
          brR = matrix(BR.R[nb, ], nrow = 1),
          cij = cij, b = b, Po = Po, ord = ord,
          Pb = matrix(BR.P[nb]),
          nbuckets = nbuckets, k = k, use_cpp = use_cpp
        )
        
        R <- rpbr$cR
        Pbr <- rpbr$pcR
        
        if (is.null(R) || length(R) == 0) next
        
        # Store all generated branches (filtering happens later in consolidation)
        valid_count <- valid_count + 1
        KR.R_list[[valid_count]] <- R
        KR.P_list[[valid_count]] <- Pbr
      }
      
      if (valid_count > 0) {
        KR.R <- do.call(rbind, KR.R_list[1:valid_count])
        KR.P <- do.call(rbind, KR.P_list[1:valid_count])
        
        # ═══════════════════════════════════════════════════════════════════
        # CRITICAL INSIGHT: Number of buckets is MONOTONICALLY INCREASING!
        # If a branch has > nbuckets now, it will NEVER decrease to nbuckets
        # Therefore: ELIMINATE all branches with > nbuckets IMMEDIATELY
        # Keep ALL branches with <= nbuckets (no limit on number)
        # ═══════════════════════════════════════════════════════════════════
        
        n_buckets_per_branch <- apply(KR.R, 1, function(x) length(unique(x)))
        valid_mask <- (n_buckets_per_branch <= nbuckets)
        
        if (any(valid_mask)) {
          # Keep only branches with <= nbuckets
          KR.R <- KR.R[valid_mask, , drop = FALSE]
          KR.P <- KR.P[valid_mask, , drop = FALSE]
          
          # Update for next iteration
          BR.R <- KR.R
          BR.P <- matrix(KR.P, ncol = 1)
        } else {
          # All branches have > nbuckets
          # Keep BR.R unchanged (don't update with invalid branches)
          if (PS) {
            message("    All ", sum(!valid_mask), " branches have > ", nbuckets, " buckets - keeping previous")
          }
          # BR.R stays as is - continue with previous valid branches
        }
      }
      
      if (PS) {
        message("  k=", k, ": ", B, " branches, ", valid_count, " valid")
      }
    }
    
    # ═══════════════════════════════════════════════════════════════════════
    # CHECK: Best k-bucket solution in this round
    # ═══════════════════════════════════════════════════════════════════════
    
    SSP <- which(BR.P == min(BR.P))
    MinP <- min(BR.P)
    
    CR_candidate <- matrix(BR.R[SSP, ], length(SSP), nobj)
    n_buckets <- apply(CR_candidate, 1, function(x) length(unique(x)))
    
    # Filter ONLY k-bucket solutions
    k_bucket_mask <- (n_buckets == nbuckets)
    
    if (PS) {
      message("  Global MinPenalty: ", round(MinP, 2), 
              " (buckets: ", paste(unique(n_buckets), collapse = ", "), ")")
    }
    
    if (any(k_bucket_mask)) {
      # Found k-bucket solutions this round
      CR_kbucket <- CR_candidate[k_bucket_mask, , drop = FALSE]
      
      # Get penalties for k-bucket solutions
      P_kbucket <- BR.P[SSP[k_bucket_mask]]
      min_kbucket_penalty <- min(P_kbucket)
      
      if (PS) {
        message("  → ", nbuckets, "-bucket penalty: ", round(min_kbucket_penalty, 2))
      }
      
      # Check if improved
      if (min_kbucket_penalty < best_kbucket_penalty) {
        # IMPROVEMENT!
        best_kbucket_penalty <- min_kbucket_penalty
        best_idx <- which(P_kbucket == min_kbucket_penalty)
        best_kbucket_solution <- CR_kbucket[best_idx, , drop = FALSE]
        no_improvement_rounds <- 0
        
        if (PS) {
          message("  ✓ IMPROVED k-bucket solution (penalty: ", 
                  round(best_kbucket_penalty, 2), ")")
        }
        
      } else {
        # No improvement
        no_improvement_rounds <- no_improvement_rounds + 1
        
        if (PS) {
          message("  No improvement (", no_improvement_rounds, " rounds)")
        }
        
        # CONVERGE if no improvement for 3 consecutive rounds
        if (no_improvement_rounds >= 3) {
          CR <- best_kbucket_solution
          if (PS) {
            message("  ✓ CONVERGED: Best ", nbuckets, "-bucket solution ",
                    "(penalty: ", round(best_kbucket_penalty, 2), ")")
          }
          WCi <- 0
          break
        }
      }
    } else {
      # No k-bucket solutions this round
      if (PS) {
        message("  → No ", nbuckets, "-bucket solutions this round ",
                "(found: ", paste(unique(n_buckets), collapse = ", "), " buckets)")
      }
      no_improvement_rounds <- no_improvement_rounds + 1
      
      # If we already have a best solution and no improvement for 5 rounds, stop
      if (!is.null(best_kbucket_solution) && no_improvement_rounds >= 5) {
        CR <- best_kbucket_solution
        if (PS) {
          message("  ✓ CONVERGED: Returning best ", nbuckets, 
                  "-bucket solution found (penalty: ", 
                  round(best_kbucket_penalty, 2), ")")
        }
        WCi <- 0
        break
      }
    }
    
    # Continue to next round
    if (WCi == 1) {
      Po <- MinP
      lambda <- lambda + 1
      nRR <- matrix(BR.R[SSP[1], ], 1, nobj)
      BR.R <- nRR
      BR.P <- 0
      a <- t(matrix(sort(BR.R)))
      ord <- t(matrix(order(BR.R)))
    }
  }
  
  # ═══════════════════════════════════════════════════════════════════════
  # FINAL CHECKS
  # ═══════════════════════════════════════════════════════════════════════
  
  if (is.null(best_kbucket_solution)) {
    # No k-bucket solution found - fallback
    warning("No ", nbuckets, "-bucket solution found during search. ",
            "Returning best approximation.")
    
    SSP <- which(BR.P == min(BR.P))
    CR <- matrix(BR.R[SSP, ], length(SSP), nobj)
    n_buckets <- apply(CR, 1, function(x) length(unique(x)))
    
    # Try to find closest to k
    closest_idx <- which.min(abs(n_buckets - nbuckets))
    CR <- matrix(CR[closest_idx, ], 1, nobj)
  } else {
    CR <- best_kbucket_solution
  }
  
  # Remove duplicates
  CR <- unique(CR)
  if (!is.matrix(CR)) CR <- matrix(CR, nrow = 1)
  
  return(CR)
}


# ══════════════════════════════════════════════════════════════════════════
# IMPROVED: Filter rankings by bucket constraint (VECTORIZED)
# ══════════════════════════════════════════════════════════════════════════

filter_by_bucket_constraint <- function(R, ord, b, nbuckets, k) {
  #' @description
  #' Improved bucket constraint check with clear logic:
  #' 
  #' ADAPTIVE CONSTRAINT STRATEGY:
  #' - Early iterations (k <= nbuckets): Allow at most nbuckets
  #' - Later iterations (k > nbuckets): Require between (nbuckets-1) and nbuckets
  #' 
  #' RATIONALE:
  #' - Early: Be permissive to explore search space
  #' - Late: Be strict to ensure final solution has exactly nbuckets
  #' 
  #' The (nbuckets-1) lower bound in late iterations allows temporary
  #' solutions that are "almost there" and will reach nbuckets when
  #' all items are assigned.
  
  if (!is.matrix(R)) {
    R <- matrix(R, nrow = 1)
  }
  
  n_candidates <- nrow(R)
  
  # Extract relevant portion of ranking (only items processed so far)
  R_subset <- R[, ord[b], drop = FALSE]
  
  # Count unique buckets for each candidate (VECTORIZED)
  n_unique_buckets <- apply(R_subset, 1, function(x) length(unique(x)))
  
  # Define valid range based on iteration
  if (k <= nbuckets) {
    # EARLY PHASE: Allow partial solutions with at most nbuckets
    valid <- (n_unique_buckets <= nbuckets)
  } else {
    # LATE PHASE: Require close to target
    # Allow (nbuckets-1) because final items may create the last bucket
    valid <- (n_unique_buckets >= (nbuckets - 1)) & (n_unique_buckets <= nbuckets)
  }
  
  if (!any(valid)) {
    return(NULL)
  }
  
  valid_idx <- which(valid)
  R_filtered <- R[valid_idx, , drop = FALSE]
  
  # Attach indices for penalty filtering
  attr(R_filtered, "valid_idx") <- valid_idx
  
  return(R_filtered)
}


# ══════════════════════════════════════════════════════════════════════════
# Phase 1 of BBconsensus (simplified for initialization)
# ══════════════════════════════════════════════════════════════════════════

BBconsensus_buckets_phase1 <- function(RR, cij, nbuckets, PS = FALSE, use_cpp = TRUE) {
  
  if (is.null(ncol(RR))) {
    RR <- matrix(RR, nrow = 1)
  }
  
  CR <- RR
  sij <- scorematrix(RR, use_cpp = use_cpp)
  Po <- sum(abs(cij)) - sum(cij * sij)
  
  a <- t(matrix(sort(RR, decreasing = TRUE)))
  ord <- t(matrix(order(RR, decreasing = TRUE)))
  R <- RR
  addpenalty <- matrix(0, length(a), 1)
  
  # Exploration of initial solution
  for (k in 2:length(a)) {
    
    b <- 1:k
    R <- ReorderingBB(R, use_cpp = use_cpp)
    KR <- t(matrix(R[ord[b]]))
    KR <- KR[-length(KR)]
    MO <- max(KR)
    MI <- min(KR)
    aa <- 1
    KO <- 1
    KR[length(KR) + 1] <- MO + 1
    R[ord[b]] <- KR
    
    candidate <- matrix(0, nrow(RR), ncol(RR))
    Pb <- matrix(0, 1, 1)
    
    while (KO == 1) {
      
      candidate <- rbind(candidate, R)
      
      if (aa == 1) {
        candidate <- matrix(candidate[-1, ], 1, ncol(candidate))
      }
      
      Sij <- scorematrix(matrix(candidate[aa, ], 1, ncol(R)), use_cpp = use_cpp)
      Pb <- rbind(Pb, sum(abs(cij)) - sum(cij * Sij))
      
      if (aa == 1) {
        Pb <- matrix(Pb[-1, ], 1, 1)
      }
      
      # IMPROVED: Use .Machine$double.xmax instead of magic number
      if (length(unique(candidate[aa, ])) != nbuckets) {
        Pb[aa] <- .Machine$double.xmax / 2  # Half of max to avoid overflow
      }
      
      if (Pb[aa] == 0) {
        CR <- R
        Po <- 0
        Pc <- 0
        break
      }
      
      Pc <- 1
      R[ord[b[length(b)]]] <- R[ord[b[length(b)]]] - 1
      
      if (MI - R[ord[b[length(b)]]] > 1) {
        KO <- 0
      }
      
      aa <- aa + 1
    }
    
    if (Pc == 0) {
      break
    }
    
    minp <- min(Pb)
    posp <- which(Pb == min(Pb))
    
    if (minp <= Po) {
      Po <- minp
      CR <- t(matrix(candidate[posp[1], ]))
      R <- CR
      addpenalty[k, 1] <- PenaltyBB2(cij, R, ord[b], use_cpp = use_cpp)
    } else {
      R <- CR
      addpenalty[k, 1] <- PenaltyBB2(cij, R, ord[b], use_cpp = use_cpp)
    }
    
    candidate <- mat.or.vec(nrow(R), ncol(R))
    Pb <- mat.or.vec(1, 1)
  }
  
  if (Pc == 0) {
    Po <- 0
    addpenalty <- 0
  } else {
    Poo <- sum(addpenalty)
  }
  
  SIJ <- scorematrix(CR, use_cpp = use_cpp)
  Po <- sum(addpenalty)
  
  return(list(cons = CR, pen = Po))
}


# ══════════════════════════════════════════════════════════════════════════
# branches_buckets: Generate branches with bucket-aware penalty filtering
# ══════════════════════════════════════════════════════════════════════════

branches_buckets <- function(brR, cij, b, Po, ord, Pb, nbuckets, k, use_cpp = TRUE) {
  
  # Generate candidate branches
  candidate <- findbranches_buckets(brR, ord, b, nbuckets, use_cpp = use_cpp)
  
  if (is.null(candidate) || nrow(candidate) == 0) {
    return(list(cR = NULL, pcR = NULL))
  }
  
  n_cand <- nrow(candidate)
  
  # Compute penalties using C++ batch if available
  if (use_cpp) {
    ord_subset <- ord[b]
    addpenalty <- PenaltyBB2_batch_impl(cij, candidate, ord_subset)
  } else {
    # R fallback
    addpenalty <- numeric(n_cand)
    for (gm in 1:n_cand) {
      addpenalty[gm] <- PenaltyBB2(cij, candidate[gm, ], ord[b], use_cpp = FALSE)
    }
  }
  
  # Compute total penalties
  Pbr <- addpenalty + as.numeric(Pb)
  
  # Filter candidates by penalty threshold
  valid_mask <- (Pbr <= Po)
  
  if (!any(valid_mask)) {
    return(list(cR = NULL, pcR = NULL))
  }
  
  # Extract valid candidates
  R <- candidate[valid_mask, , drop = FALSE]
  Pbr_valid <- Pbr[valid_mask]
  
  # Ensure matrix format
  if (nrow(R) == 1) {
    R <- matrix(R, 1, ncol(candidate))
  }
  
  return(list(cR = R, pcR = matrix(Pbr_valid, ncol = 1)))
}


# ══════════════════════════════════════════════════════════════════════════
# findbranches_buckets: Generate branches with native bucket constraint
# ══════════════════════════════════════════════════════════════════════════

findbranches_buckets <- function(R, ord, b, nbuckets, use_cpp = TRUE) {
  
  N <- ncol(R)
  k <- length(b)
  
  # ═══════════════════════════════════════════════════════════════
  # STEP 1: Analizza situazione corrente
  # ═══════════════════════════════════════════════════════════════
  
  R_partial <- R[ord[b]]  # Ranking parziale
  current_buckets <- length(unique(R_partial[-length(R_partial)]))  # Escludi ultimo
  remaining_items <- N - k
  
  # ═══════════════════════════════════════════════════════════════
  # STEP 2: Calcola quanti buckets DEVI ancora creare
  # ═══════════════════════════════════════════════════════════════
  
  buckets_needed <- nbuckets - current_buckets
  
  # ═══════════════════════════════════════════════════════════════
  # STEP 3: Verifica fattibilità
  # ═══════════════════════════════════════════════════════════════
  
  if (buckets_needed > remaining_items) {
    # Impossibile: troppi buckets da creare con pochi items
    return(NULL)
  }
  
  if (buckets_needed < 0) {
    # Impossibile: già troppi buckets
    return(NULL)
  }
  
  # ═══════════════════════════════════════════════════════════════
  # STEP 4: Genera candidati VALIDI (FIXED seq() logic)
  # ═══════════════════════════════════════════════════════════════
  
  last_rank <- R_partial[length(R_partial)]
  min_rank <- min(R_partial[-length(R_partial)])
  
  if (buckets_needed == 0) {
    # SOLO tie con bucket esistenti
    valid_ranks <- unique(R_partial[-length(R_partial)])
    
  } else if (buckets_needed == remaining_items) {
    # DEVI creare un nuovo bucket per OGNI remaining item
    start_val <- last_rank - 1
    end_val <- min_rank - 1
    
    # FIX: Gestisci caso start < end
    if (start_val >= end_val) {
      valid_ranks <- seq(start_val, end_val, by = -1)
    } else {
      valid_ranks <- start_val  # Almeno un nuovo bucket
    }
    
  } else {
    # Flessibilità: tie O nuovo bucket
    start_val <- last_rank - 1
    end_val <- min_rank - 1
    
    # FIX: Gestisci caso start < end
    if (start_val >= end_val) {
      new_buckets <- seq(start_val, end_val, by = -1)
    } else {
      new_buckets <- start_val  # Almeno un nuovo bucket
    }
    
    valid_ranks <- c(last_rank, new_buckets)
  }
  
  # ═══════════════════════════════════════════════════════════════
  # STEP 5: Genera matrice candidati
  # ═══════════════════════════════════════════════════════════════
  
  n_candidates <- length(valid_ranks)
  candidates <- matrix(0, n_candidates, N)
  
  for (i in 1:n_candidates) {
    candidates[i, ] <- R
    candidates[i, ord[k]] <- valid_ranks[i]
  }
  
  return(candidates)
}

#-------------------------------------