
B <- matrix(c(50, 67, 90, 98, 120,
             55, 71, 93, 102, 129,
             65, 76, 95, 105, 134,
             50, 80, 102, 130, 138,
             60, 82, 97, 135, 151,
             65, 89, 106, 137, 153,
             75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
rownames(B) <- c("G1","G2","G3","G4","G5","G6","G7")
colnames(B) <- c("E1","E2","E3","E4","E5")
 
B2 = B
B2[1,1] = B2[2,1] = NA

m4 = RnipalsPca(scale(B2), nPcs=5, threshold=1e-9)

# P loadings, do NOT match ade4, plsdepot, mixOmics
m4$loadings
round(t(m4$loadings) %*% m4$loadings,3)

# T scores
m4$scores

# ----- gge -----

m5 = gge::nipals(scale(B2), maxcomp=5)
m5$eval
sum(m5$eval) # 27.73141

# P loadings, same as model m4 from pcaMethods

m5$rotation
# T scores
m5$scores

RnipalsPca <- function(Matrix, nPcs=2, 
                       varLimit=1,
                       maxSteps=5000, 
                       threshold=1e-6, verbose=interactive(),
                       ##' PCA by non-linear iterative partial least squares
##'
##' Can be used for computing PCA on a numeric matrix using either the
##' NIPALS algorithm which is an iterative approach for estimating the
##' principal components extracting them one at a time. NIPALS can
##' handle a small amount of missing values. It is not recommended to
##' use this function directely but rather to use the pca() wrapper
##' function. There is a C++ implementation given as \code{nipalsPca}
##' which is faster.
##' @title NIPALS PCA implemented in R
##' @param Matrix Pre-processed (centered, scaled) numerical matrix
##' samples in rows and variables as columns. 
##' @param nPcs Number of components that should be extracted.
##' @param varLimit Optionally the ratio of variance that should be
##' explained. \code{nPcs} is ignored if varLimit < 1
##' @param maxSteps Defines how many iterations can be done before
##' algorithm should abort (happens almost exclusively when there were
##' some wrong in the input data).
##' @param threshold The limit condition for judging if the algorithm
##' has converged or not, specifically if a new iteration is done if
##' \eqn{(T_{old} - T)^T(T_{old} - T) > \code{limit}}.
##' @param verbose Show simple progress information.
##' @param ... Only used for passing through arguments.
##' @return A \code{pcaRes} object.
##' @references Wold, H. (1966) Estimation of principal components and
##' related models by iterative least squares. In Multivariate
##' Analysis (Ed., P.R. Krishnaiah), Academic Press, NY, 391-420.
##' @author Henning Redestig
##' @seealso \code{prcomp}, \code{princomp}, \code{pca}
##' @examples
##' data(metaboliteData)
##' mat <- prep(t(metaboliteData))
##' ## c++ version is faster
##' system.time(pc <- RnipalsPca(mat, method="rnipals", nPcs=2))
##' system.time(pc <- nipalsPca(mat, nPcs=2))
##' ## better use pca()
##' pc <- pca(t(metaboliteData), method="rnipals", nPcs=2)
##' \dontshow{stopifnot(sum((fitted(pc) - t(metaboliteData))^2, na.rm=TRUE) < 200)}
##' @keywords multivariate
##' @export
RnipalsPca <- function(Matrix, nPcs=2, 
                       varLimit=1,
                       maxSteps=5000, 
                       threshold=1e-6, verbose=interactive(), ...) {

  nVar <- ncol(Matrix)

  ##Find a good? starting column -- better way?
  startingColumn <- 1

  ## sum(c(NA, NA), na.rm=TRUE) is 0, but we want NA
  sum.na <- function(x){ ifelse(all(is.na(x)), NA, sum(x, na.rm=TRUE))}

  TotalSS <- sum(Matrix*Matrix, na.rm=TRUE)

  ph <- rep(0, nVar)
  R2cum <- rep(NA, nPcs)
  scores <- NULL
  loadings <- NULL
  anotherPc <- TRUE
  l <- 1
  
  while(anotherPc) {
    count <- 0                 #number of iterations done
    th <- Matrix[,startingColumn]   #first column is starting vector for th
    continue <- TRUE
    if(verbose) cat(paste("Calculating PC", l, ": ", sep=""))
    
    while(continue) {
      count <- count+1
      ph <- rep(0, nVar)

      ##Calculate loadings through LS regression
      ##Note: Matrix*th is column-wise multiplication
      tsize <- sum(th * th, na.rm=TRUE)
      ph <- apply(Matrix * (th / tsize), 2, sum.na)
      ##normalize ph based on the available values.
      psize <- sum(ph*ph, na.rm=TRUE)
      ph <- ph / sqrt(psize)
      
      ##Calculate scores through LS regression
      ##Trick: To get row-wise multiplication, use t(Matrix)*ph, then
      ##be sure to use apply(,2,) and NOT apply(,1,)!
      th.old <- th
      th <- apply(t(Matrix) * ph, 2, sum.na)
      
      ##Round up by calculating if convergence condition is met and
      ##checking if it seems to be an neverending loop.
      if (count > maxSteps) {
        stop("Too many iterations, quitting")
      }
      if (t(na.omit(th.old - th)) %*% (na.omit(th.old - th)) <= threshold) {
        continue = FALSE
      }
      if (verbose)cat("*")
    }
    if (verbose) cat(" Done\n")
    Matrix <- Matrix - (th %*% t(ph))
    scores <- cbind(scores, th)
    loadings <- cbind(loadings, ph)
    
    ##cumulative proportion of variance
    R2cum[l] <- 1 - (sum(Matrix*Matrix,na.rm=TRUE) / TotalSS)
    l <- l + 1
    if((!abs(varLimit - 1) < 1e-4 & R2cum[l - 1] >= varLimit) | l > nPcs) {
      anotherPc <- FALSE
      nPcs <- l - 1
    }
  }

  #r <- new("pcaRes")
  return(list(
    scores = scores,
    loadings = loadings,
    R2cum = R2cum,
    varLimit = varLimit,
    method = "rnipals"))

}
