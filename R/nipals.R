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
#' The R2 values that are reported are marginal, not cumulative.
#' 
#' @param x Numerical matrix for which to find principal compontents. 
#' Missing values are allowed.
#' 
#' @param ncomp Maximum number of principal components to extract from x.
#'
#' @param center If TRUE, subtract the mean from each column of x before PCA.
#'
#' @param scale if TRUE, divide the standard deviation from each column of x before PCA.
#' 
#' @param maxiter Maximum number of NIPALS iterations for each
#' principal component.
#'
#' @param tol Default 1e-6 tolerance for testing convergence of the NIPALS
#' iterations for each principal component.
#'
#' @param startcol Determine the starting column of x for the iterations
#' of each principal component.
#' If 0, use the column of x that has maximum absolute sum.
#' If a number, use that column of x.
#' If a function, apply the function to each column of x and choose the column
#' with the maximum value of the function.
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
#' \code{fitted}, \code{ncomp}, \code{R2}, \code{iter}, \code{center}, 
#' \code{scale}.
#' 
#' @references
#' Wold, H. (1966) Estimation of principal components and
#' related models by iterative least squares. In Multivariate
#' Analysis (Ed., P.R. Krishnaiah), Academic Press, NY, 391-420.
#'
#' Andrecut, Mircea (2009).
#' Parallel GPU implementation of iterative PCA algorithms.
#' Journal of Computational Biology, 16, 1593-1599.
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
#' # Two ways to make a biplot
#'
#' # method 1
#' biplot(p2$scores, p2$loadings)
#'
#' # method 2
#' class(p2) <- "princomp"
#' p2$sdev <- sqrt(p2$eig)
#' biplot(p2, scale=0)
#' 
#' @author Kevin Wright
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
  nvar <- ncol(x)
  nobs <- nrow(x)
  x.orig <- x # Save x for row/col names

  # Check for a column or row with all NAs
  col.na.count <- apply(x, 2, function(x) sum(!is.na(x)))
  if(any(col.na.count==0)) stop("At least one column is all NAs")
  row.na.count <- apply(x, 1, function(x) sum(!is.na(x)))
  if(any(row.na.count==0)) stop("At least one row is all NAs")
  
  # center / scale
  if(center) {
    cmeans <- colMeans(x, na.rm=TRUE)
    x <- sweep(x, 2, cmeans, "-")
  } else cmeans <- NA
  if(scale) {
    csds <- apply(x, 2, sd, na.rm=TRUE)
    x <- sweep(x, 2, csds, "/")
  } else csds <- NA
  
  TotalSS <- sum(x*x, na.rm=TRUE)
  
  # initialize outputs
  PPp = matrix(0, nrow=nvar, ncol=nvar)
  TTp = matrix(0, nrow=nobs, ncol=nobs)
  eig <- rep(NA, length=ncomp)
  R2cum <- rep(NA, length=ncomp)
  loadings <- matrix(nrow=nvar, ncol=ncomp)
  scores <- matrix(nrow=nobs, ncol=ncomp)
  iter <- rep(NA, length=ncomp)

  # position of NAs
  x.miss <- is.na(x)
  has.na <- any(x.miss)
  if(force.na) has.na <- TRUE

  # Calculate PC h
  for(h in 1:ncomp) {
    
    # start with a column of x, call it th
    if(is.function(startcol)){
      scol <- which.max(apply(x, 2, startcol))
    } else if(startcol==0L) {
      # use the column with maximum absolute sum
      scol <- which.max(apply(x, 2, function(x) sum(abs(x), na.rm=TRUE)) )
      #scol <- which.max(apply(x, 2, var, na.rm = TRUE))
    } else {
      scol <- startcol
    }
    if(verbose >= 1) cat("PC ", h, " starting column: ", scol, sep="")

    # replace NA values with 0 so those elements don't contribute
    # to dot-products, etc
    if(has.na){
      x0 <- x
      x0[x.miss] <- 0
      th <- x0[, scol]
    } else {
      th <- x[, scol]
    }
    
    pciter <- 1 # reset iteration counter for each PC
    continue <- TRUE
    while(continue) {
  
      # loadings p = X't/t't  
      if(has.na){
        # caution: t't is NOT the same for each column of X't, but is
        # the sum of the squared elements of t for which that column
        # of X is not missing data
        T2 <- matrix(th*th, nrow=nobs, ncol=nvar)
        T2[x.miss] <- 0
        # it sometimes happen that colSums(T2) has a 0. should we check
        # for that or just let it fail?
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
        P2 <- matrix(ph*ph, nrow=nvar, ncol=nobs)
        P2[t(x.miss)] <- 0
        th = x0 %*% ph / colSums(P2)        
      } else {
        th = x %*% ph / sum(ph*ph)
      }

      # Gram Schmidt orthogonalization # t = t - (Th)(Th)' t
      if(gramschmidt && h>1) {
        th <- th - TTp %*% th
      }
      
      # check convergence of th
      if( sum((th-th.old)^2, na.rm=TRUE)<tol ) continue=FALSE
      
      pciter <- pciter + 1
      if (pciter == maxiter) {
        continue <- FALSE
        warning("Stopping after ", maxiter, " iterations for PC ", h,".\n")
      }
      
      if (verbose >= 1) cat(".")
    } # iterations for PC h
    if (verbose >= 1) cat("\n")
    
    # deflate/remove variation from x explained by PC h, x=x-tp' 
    x <- x - (th %*% t(ph))
    loadings[,h] <- ph
    scores[,h] <- th
    eig[h]  <- sum(th*th, na.rm=TRUE)
    iter[h] <- pciter
    
    # Update (Ph)(Ph)' and (Th)(Th)' for next PC
    if(gramschmidt) {
      PPp = PPp + tcrossprod(ph)          # PP' = PP' + (ph)(ph)'
      TTp = TTp + tcrossprod(th) / eig[h] # TT' = TT' = (th)(th)'
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
    # must use diag( nrow=length(eig)) because diag(3.3) is a 3x3 identity
    # xhat <- tcrossprod( tcrossprod(scores,diag(eig, nrow=length(eig))), loadings)
    xhat <- tcrossprod( sweep( scores, 2, eig, "*") , loadings)
    if(scale) xhat <- sweep(xhat, 2, csds, "*")
    if(center) xhat <- sweep(xhat, 2, cmeans, "+")
    rownames(xhat) <- rownames(x.orig)
    colnames(xhat) <- colnames(x.orig)
  } else {
    xhat <- NULL
  }
  
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste("PC", 1:ncol(scores), sep="")
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- paste("PC", 1:ncol(loadings), sep="")
  
  out <- list(eig=eig,
              scores=scores, 
              loadings=loadings,
              fitted=xhat,
              ncomp=ncomp,
              R2=R2,
              iter=iter, 
              center=cmeans, scale=csds)
  return(out)
}

## ---------------------------------------------------------------------------


#' Average angular distance between two rotation matrices
#'
#' For matrices A and B, calculate the angle between the column
#' vectors of A and the corresponding column vectors of B.
#' Then average the angles.
#'
#' The results of the singular value decomposition X=USV' are
#' unique, but only up to a change of sign for columns of U,
#' which indicates that the axis is flipped.
#' 
#' @param A Matrix
#' 
#' @param B Matrix
#' 
#' @return
#' A single floating point number, in radians.
#' 
#' @author Kevin Wright
#' 
#' @examples
#' # Example from https://math.stackexchange.com/questions/2113634/
#' rot1 <- matrix(c(-0.956395958, -0.292073218, 0.000084963,
#'                  0.292073230, -0.956395931, 0.000227268,
#'                  0.000014880, 0.000242173, 0.999999971),
#'                ncol=3, byrow=TRUE)
#' rot2 <- matrix(c(-0.956227882, -0.292623029, -0.000021887,
#'                  0.292623030, -0.956227882, -0.000024473,
#'                  -0.000013768, -0.000029806, 0.999999999),
#'                ncol=3, byrow=TRUE)
#' avg_angular_distance(rot1, rot2) # .0004950387
#' 
#' @references
#' None
#' 
#' @export 
avg_angular_distance <- function(A, B){
  dotprod <- colSums(A*B) / sqrt( colSums(A*A) * colSums(B*B) )
  
  # Note: Use abs() so that vectors which point in nearly opposite
  # directions ( dotprod= -1 ) can be considered as if pointing
  # in the same direction.
  # Note: we need to make sure elements of dotprod are not bigger
  # than 1.0, which sometimes happens due to floating point error.
  # theta[i] is the angle between A[,i] and B[,i]
  #theta <- acos( abs(dotprod) - .Machine$double.eps )
  theta <- acos( pmin( abs(dotprod), 1 ) )
  # average radian angle between pair-wise columns of A and B
  avg_angle <- mean(theta)
  return(avg_angle)
}

