#' Compute Randles Interdirections
#'
#' Computes the matrix of Randles interdirections for multivariate data
#' in dimension \eqn{d \ge 3}.
#'
#' For observations \eqn{X_1,\dots,X_n \in \mathbb{R}^d}, the
#' \eqn{(i,j)} entry is the number of hyperplanes through the origin
#' determined by \eqn{d-1} observations other than \eqn{X_i} and
#' \eqn{X_j} that separate \eqn{X_i} and \eqn{X_j}.
#'
#' @param X A numeric matrix with observations in rows and variables
#'   in columns. The dimension must satisfy \code{ncol(X) >= 3}.
#' @param progress Logical. If \code{TRUE}, display progress information
#'   for dimensions greater than three.
#' @param tol Numerical tolerance used for rank and degeneracy checks.
#'
#' @return A symmetric \eqn{n \times n} matrix of Randles
#'   interdirection counts, with zeros on the diagonal.
#'
#' @details
#' The implementation assumes the data are in general position.
#'
#' In dimension three, the algorithm uses a gnomonic projection and a
#' sweep of the associated planar line arrangement.
#'
#' For dimensions greater than three, the problem is reduced to
#' three-dimensional problems by quotienting out spans of
#' \eqn{d-3} observations.
#'
#' With the current three-dimensional implementation, which sorts the
#' arrangement crossings, the running time for fixed dimension \eqn{d}
#' is \eqn{O(n^{d-1}\log n)}.
#'
#' @examples
#' set.seed(123)
#' X <- matrix(rnorm(40), ncol = 4)
#' C <- randles_interdirections(X)
#' C
#'
#' @export
randles_interdirections <- function(X, progress = FALSE, tol = 1e-10) {

  X <- as.matrix(X)

  if (!is.numeric(X)) {
    stop("X must be a numeric matrix.")
  }

  if (ncol(X) < 3L) {
    stop("The dimension must satisfy ncol(X) >= 3.")
  }

  if (any(!is.finite(X))) {
    stop("X must contain only finite values.")
  }

  randles_general_rcpp(
    X,
    progress = progress,
    tol = tol
  )
}
