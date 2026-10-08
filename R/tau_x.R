#' TauX (tau exstension) rank correlation coefficient
#'
#' Tau exstension is a new rank correlation coefficient defined by Emond and Mason (2002)
#'
#' @param X a M by N data matrix, in which there are N judges and M objects to be judged. Each row is a ranking of the objects which are represented by the columns. If there is only X as input, the output is a square matrix containing the Tau_X rcc.
#' @param Y A row vector, or a n by M data matrix in which there are n judges and the same M objects as X to be judged.
#'
#' @return Tau_x rank correlation coefficient
#'
#' @details With complete rankings TauX is a linear function of the Kemeny distance. The rankings can also be partial (missing objects coded NA): as in Emond and Mason (2002), the score of a pair of objects with a missing object is 0, so that the pair does not contribute to the coefficient.
#'
#' @references Emond, E. J., and Mason, D. W. (2002). A new rank correlation coefficient with application to the consensus ranking problem. Journal of Multi-Criteria Decision Analysis, 11(1), 17-28.
#' 
#' @author Antonio D'Ambrosio \email{antdambr@unina.it}
#' 
#' @seealso \code{\link{kemenyd}} Kemeny distance
#' @seealso \code{\link{iw_tau_x}} item-weighted tau_x rank correlation coefficient
#' 
#' @keywords TauX rank correlation coefficient
#' 
#' @examples
#' data(BU)
#' RD<-BU[,1:3]
#' Tau<-tau_x(RD)
#' Tau1_3<-tau_x(RD[1,],RD[3,])
#' 
#' @export 


tau_x <- function(X,Y=NULL) {
  
  if(is(X,"numeric") & !is(X,"matrix")){
    X<-matrix(X,ncol=length(X))
  }
  
  
  N<-ncol(X)
  maxd<-N*(N-1)
  if (!is(Y,"NULL") && is(Y,"numeric") & !is(Y,"matrix")){
    Y<-matrix(Y,ncol=length(Y))
  }
  
  if (!anyNA(X) && (is(Y,"NULL") || !anyNA(Y))) {
    # complete rankings: TauX from the Kemeny distance
    if (is(Y,"NULL")){
      d<-kemenyd(X)
      tau <- 1-(2*d/maxd)
    } else {
      d<-kemenyd(X,Y)
      tau <- 1-(2*d/maxd)
    }
    return(tau)
  }
  
  # partial rankings (NA): TauX as defined by Emond and Mason (2002), from the
  # score matrices, in which the pairs with a missing object are coded 0. The
  # result has the same form as for complete rankings
  scorevec <- function(Z) {
    t(apply(Z, 1, function(r) as.vector(scorematrix(r, use_cpp = FALSE))))
  }
  SX <- scorevec(X)
  if (is(Y,"NULL")){
    tau <- 1-(2*kemenyd(matrix(1, nrow(X), N, dimnames = dimnames(X)))/maxd)
    TT <- tcrossprod(SX) / maxd
    tau[] <- TT[lower.tri(TT)]
  } else {
    tau <- 1-(2*kemenyd(matrix(1, nrow(X), N, dimnames = dimnames(X)),
                        matrix(1, nrow(Y), N, dimnames = dimnames(Y)))/maxd)
    tau[] <- tcrossprod(SX, scorevec(Y)) / maxd
  }
  tau 
}

#' @rdname tau_x
#' @export
Tau_X <- tau_x

