# nipals.R

if(FALSE){

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

m4 <- nipals(B2, ncomp=5)
m4$eig

}

#' Principal component analysis by NIPALS, non-linear iterative partial least squares
#'
#' Used for finding principal components of a numeric matrix.
#' Missing values in the matrix are allowed.
#' Principal Components are extracted one a time.  
#' The algorithm computes x = TP', where T is the 'scores' matrix and P is
#' the 'loadings' matrix.
#' 
#' @param x Numerical matrix for which to find principal compontents. 
#' Missing values are allowed.
#' 
#' @param ncomp Maximum number of principal components to extract from x.
#'
#' @param center If TRUE, subtract the mean from each column of x.
#'
#' @param scale if TRUE, divide the standard deviation from each column of x.
#' 
#' @param maxiter Maximum number of NIPALS iterations for each
#' principal component.
#'
#' @param tol Default 1e-9 tolerance for testing convergence of the NIPALS
#' iterations for each principal component.
#'
#' @param startcol If 0, start the iterations for each principal component
#' with the column of x that has maximum variation.
#' Otherwise, start with the spcified column number.
#'
#' @param fitted Default FALSE. If TRUE, return the fitted (reconstructed) value of x.
#'
#' @param force.na Default FALSE. If TRUE, force the function to use the
#' method for missing values, even if there are no missing values in x.
#' 
#' @param gramschmidt Default TRUE. If TRUE, perform Gram-Schmidt 
#' orthogonalization at each iteration.
#' 
#' @param verbose Default FALSE. Use TRUE or 1 to show some diagnostics.
#' 
#' @return A list with components \code{eig}, \code{scores}, \code{loadings}, 
#' \code{ncomp}, \code{R2}, \code{xhat}.
#' 
#' @references
#' Wold, H. (1966) Estimation of principal components and
#' related models by iterative least squares. In Multivariate
#' Analysis (Ed., P.R. Krishnaiah), Academic Press, NY, 391-420.
#' 
#' @examples 
#' B <- matrix(c(50, 67, 90, 98, 120,
#'               55, 71, 93, 102, 129,
#'               65, 76, 95, 105, 134,
#'               50, 80, 102, 130, 138,
#'               60, 82, 97, 135, 151,
#'               65, 89, 106, 137, 153,
#'               75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
#' rownames(B) <- c("G1","G2","G3","G4","G5","G6","G7")
#' colnames(B) <- c("E1","E2","E3","E4","E5")
#' dim(B) # 7 x 5
#' p1 <- nipals(B)
#' dim(p1$scores) # 7 x 5
#' dim(p1$loadings) # 5 x 5
#' 
#' B2 = B
#' B2[1,1] = B2[2,2] = NA
#' p2 = nipals(B2, fitted=TRUE)
#' 
#' @author Kevin Wright.
#' 
#' @importFrom stats sd var
#' @export
nipals <- function(x,
                   ncomp=min(nrow(x), ncol(x)),
                   center=TRUE, scale=TRUE,
                   maxiter=500,
                   tol=1e-6,
                   startcol=0,
                   fitted=FALSE,
                   force.na=FALSE,
                   gramschmidt=TRUE,
                   verbose=FALSE) {

  x <- as.matrix(x) # in case it is a data.frame
  nc <- ncol(x)
  nr <- nrow(x)
  x.orig <- x # Save x for replacing missing values

  # Check for a column or row with all NAs
  col.count <- apply(x, 2, function(x) sum(!is.na(x)))
  if(any(col.count==0)) stop("At least one column is all NAs")
  row.count <- apply(x, 1, function(x) sum(!is.na(x)))
  if(any(row.count==0)) stop("At least one row is all NAs")
  
  # center / scale
  if(center) {
    cmeans <- colMeans(x, na.rm=TRUE)
    x <- sweep(x, 2, cmeans, "-")
  }
  if(scale) {
    csds <- apply(x, 2, sd, na.rm=TRUE)
    x <- sweep(x, 2, csds, "/")
  }
  
  TotalSS <- sum(x*x, na.rm=TRUE)
  
  # initialize outputs
  PPp = matrix(0, nrow=nc, ncol=nc)
  TTp = matrix(0, nrow=nr, ncol=nr)
  eig <- rep(NA, length=ncomp)
  R2cum <- rep(NA, length=ncomp)
  loadings <- matrix(nrow=nc, ncol=ncomp)
  scores <- matrix(nrow=nr, ncol=ncomp)

  # position of NAs
  x.miss <- is.na(x)
  has.na <- any(x.miss)
  if(force.na) has.na <- TRUE # 

  # Calculate PC h
  for(h in 1:ncomp) {
    
    # take a column of x, call it th
    # if startcol=0, use the column with maximum variance
    if(startcol==0L) {
      scol <- which.max(apply(x,2,var,na.rm = TRUE))
    } else {
      scol <- startcol
    }
    if(verbose >= 1) cat("PC ", h, " starting column: ", scol, sep="")
    th <- x[, scol]
    
    # replace NA values with 0 so those elements don't contribute
    # to dot-products, etc
    if(has.na){
      x0 <- x
      x0[x.miss] <- 0
      th <- x0[, scol]
    } else {
      th <- x[, scol]
    }
    
    iter <- 1 # reset iteration counter for each PC
    continue <- TRUE
    while(continue) {
  
      # loadings p = X't/t't  
      if(has.na){
        # caution: t't is NOT the same for each column of X't, but is
        # the sum of the squared elements of t for which that column
        # of X is not missing data
        T2 <- matrix(th*th, nrow=nr, ncol=nc)
        T2[x.miss] <- 0
        ph <- crossprod(x0,th) / colSums(T2)
      } else {
        ph = crossprod(x,th) / sum(th*th)
      }
      
      # Gram Schmidt orthogonalization p = p - PhPh'p
      if(gramschmidt && h>1) {
        ph <- ph - PPp %*% ph
      }
      # normalize to unit length p = p / p'p
      ph <- ph / sqrt(sum(ph*ph, na.rm=TRUE))
      
      # scores t = Xp/p'p
      th.old <- th
      if (has.na) {
        # square the elements of the vector ph, put into columns of P2,
        # extract the non-missing (in each column of X), and sum
        P2 <- matrix(ph*ph, nrow=nc, ncol=nr)
        P2[t(x.miss)] <- 0
        th = x0 %*% ph / colSums(P2)        
      } else {
        th = x %*% ph / sum(ph*ph)
      }

      # Gram Schmidt # t = t - T_h T_h' t
      if(gramschmidt && h>1) {
        th <- th - TTp %*% th
      }
      
      # check convergence
      if (iter > maxiter) {
        continue <- FALSE
        warning("Stopping after ", maxiter, " iterations for PC ", h,".\n")
      }
      if( sum((th-th.old)^2, na.rm=TRUE)<tol ) continue=FALSE
      
      if (verbose >= 1) cat(".")
      iter <- iter+1
    } # iterations for PC h
    if (verbose >= 1) cat("\n")
    
    # deflate/remove variation from x explained by PC h, x=x-tp' 
    x <- x - (th %*% t(ph))
    loadings[,h] <- ph
    scores[,h] <- th
    eig[h] = sum(th*th, na.rm=TRUE)
    
    # Update (Ph)(Ph)' and (Th)(Th)' for next PC
    if(gramschmidt) { 
      PPp = PPp + tcrossprod(ph)
      TTp = TTp + tcrossprod(th) / eig[h]
    }
    
    # Cumulative proportion of variance explained
    R2cum[h] <- 1 - (sum(x*x,na.rm=TRUE) / TotalSS)
    
  } # Done finding PCs
  
  # un-cumulate R2
  R2 <- c(R2cum[1], diff(R2cum))
  
  # sweep out eigenvalues from scores
  eig = sqrt(eig)
  scores = sweep(scores, 2, eig, "/")

  if(fitted) {
    # re-construction of x using ncomp principal components
    xhat <- tcrossprod( tcrossprod(scores,diag(eig)), loadings)
    if(scale) xhat <- sweep(xhat, 2, csds, "*")
    if(center) xhat <- sweep(xhat, 2, cmeans, "+")
    rownames(xhat) <- rownames(x.orig)
    colnames(xhat) <- colnames(x.orig)
  } else { xhat <- NULL }
  
  # output
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste("PC", 1:ncol(scores), sep="")
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- paste("PC", 1:ncol(loadings), sep="")
  
  out <- list(eig=eig,
              scores=scores, 
              loadings=loadings,
              fitted=xhat,
              ncomp=ncomp,
              R2=R2)
  return(out)
}
