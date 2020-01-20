# empca.R

#' Principal component analysis by weighted EMPCA, expectation maximization 
#' principal component-analysis
#' 
#' Used for finding principal components of a numeric matrix.
#' Missing values in the matrix are allowed.
#' Weights for each element of the matrix are allowed.
#' Principal Components are extracted one a time.  
#' The algorithm computes x = TP', where T is the 'scores' matrix and P is
#' the 'loadings' matrix.
#' 
#' @param x Numerical matrix for which to find principal components. 
#' Missing values are allowed.
#' 
#' @param w Numerical matrix of weights.
#' 
#' @param ncomp Maximum number of principal components to extract from x.
#'
#' @param center If TRUE, subtract the mean from each column of x before PCA.
#'
#' @param scale if TRUE, divide the standard deviation from each column of x before PCA.
#' 
#' @param maxiter Maximum number of EM iterations for each
#' principal component.
#'
#' @param tol Default 1e-6 tolerance for testing convergence of the EM
#' iterations for each principal component.
#'
#' @param seed Random seed to use when initializing the random rotation matrix.
#' 
#' @param fitted Default FALSE. If TRUE, return the fitted (reconstructed) value of x.
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
#' Stephen Bailey (2012). 
#' Principal Component Analysis with Noisy and/or Missing Data.
#' Publications of the Astronomical Society of the Pacific.
#' http://doi.org/10.1086/668105
#' 
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
#' p1 <- empca(B)
#' dim(p1$scores) # 7 x 5
#' dim(p1$loadings) # 5 x 5
#' 
#' B2 = B
#' B2[1,1] = B2[2,2] = NA
#' p2 = empca(B2, fitted=TRUE)
#' 
#' @author Kevin Wright
#' 
#' @importFrom stats rnorm lm.wfit sd var
#' @export
empca <- function(x, w,
                  ncomp=min(nrow(x), ncol(x)),
                  center=TRUE, scale=TRUE,
                  maxiter=100,
                  tol=1e-6,
                  seed=NULL,
                  fitted=FALSE,
                  gramschmidt=TRUE,
                  verbose=FALSE) {

  x <- as.matrix(x) # in case it is a data.frame
  nvar <- ncol(x)
  nobs <- nrow(x)
  x.orig <- x # Save x for replacing missing values
  
  if(missing(w)) {
    w = !is.na(x) # force missing values to have weight 0 
  }
  
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
  
  # position of NAs
  x.miss <- is.na(x)
  has.na <- any(x.miss)
  # Change missing values to 0 with 0 weight
  x[x.miss] <- 0
  w[x.miss] <- 0
  
  
  # intialize C,P
  C <- matrix(NA, nrow=nobs, ncol=ncomp)
  
  if(missing(seed)) {
    # just use identity matrix for starting rotation
    P <- diag(max(ncomp, nvar))
  } else {
    # P matrix with random orthonormal columns, nvar*ncomp
    set.seed(seed)
    ranmat <- matrix(rnorm(nobs*ncomp), nrow=nobs, ncol=ncomp)
    #qrmod <- qr(ranmat)
    #di <- diag(qr.R(qrmod))
    # di <- di/abs(di) # di is vector of +1 and -1. Why needed???
    #qmat <- qr.Q(qrmod)
    #P <- qmat %*% diag( di/abs(di) ) %*% qmat
    P <- qr.Q(qr(ranmat))
  }
  P <- P[1:nvar, 1:ncomp]
  # zapsmall( t(P) %*% P ) # P is orthonormal in columns

  for(emiter in 1:maxiter) {
    if(verbose >  1) cat(emiter,"/", maxiter, "----------------\n")
    P0 <- P # previous P to check convergence
    
    # E step
    # Find C such that X=CP
    # Find C[i,k] such that X[i,] = Sum_k: C[i,k] P[,k]
    C[] <- 0 # reset to 0  
    for(i in 1:nobs) {  
      # C[i,] <- t( solve_weighted(P, x[i,], w[i,]) )
      C[i,] <- t(lm.wfit(P, x[i,], w[i,])$coef)
      # x - replace(C, is.na(C), 0) %*% P
    }
    #browser()
    # replace NA by 0
    C[is.na(C)] <- 0
    if(verbose > 1) {
      cat("C\n")
      print(signif(C,2))
    }
    if(verbose > 1) {
      cat("x-CP\n")
      print(round(x- C %*% P,2))
    }
    
    # M step: Calculate eigvectors in P
    x1 = x # reset for each EM iteration
    for( h in 1:ncomp ) {
      # Bailey eqn 21.
      # P[j,h] = sum_j (x[,j] * w[,j] * C[,h]) / sum_j (w[,j]*C[,h]*C[,h])
      # P[ ,h] = colSums(x * w *C[,h]) / colSums(w * C[,h] * C[,h])
      wC = w*C[,h] # column C[,h] times each column of w
      P[,h] <- colSums(x1 * wC) / colSums(wC * C[,h])
      # If a column of C is entirely NA, then P[,h] is NA. Should we just
      # set it 0 instead?
      if(all(is.na(P[,h]))) P[,h] <- 0
      # not necessary to unitize columns P[,h]
      # x1 <- x1 - P[,h] %*% t(C[,h]) # deflate, Bailey eqn 22.
      x1 <- x1 - outer(C[,h],P[,h])
      x1[x.miss] <- 0 # 
    }
    
    if(gramschmidt) {
      # Re-normalize and re-orthogonalize columns of P via GramSchmidt
      # https://stackoverflow.com/questions/15584221/gram-schmidt-with-r
      # QR is not unique...flip the sign of any column of Q and row of R
      # When we do the Gram-Schmidt step on C and (separately) on P, we
      # could end up flipping C[,1] differently than P[,1]. To correct this,
      # for each i where R[i,i] < 0, flip sign of of R[i,] and Q[,i] so that
      # all R[i,i] are positive.
      # https://math.stackexchange.com/questions/2237262/
      
      # extract eigenvalues before C is orthonormalized
      eig <- sqrt(colSums(C^2))
      
      #browser()
      # perform gram-schmidt rotation on C
      qrc <- qr(C)
      C <- qr.Q(qrc) # C is now orthonormal with unit-length columns
      # flip sign of columns of C where diagonal elements of R are negative
      # do NOT use these 2 lines; fails if any di==0
      #di <- diag(qr.R(qrc))
      #di <- di/abs(di) # vector of +1 and -1
      di <- ifelse( diag(qr.R(qrc))<0, -1, 1)
      C <- C * rep(di, each=nrow(C)) # flip sign of columns as needed
      # zapsmall( crossprod(C)) # prove C'C orthonormal
      
      # do we really need to do this, or could we get information from the
      # gram-schmidt rotation of C ???
      
      # gram-schmidt on P
      qrp <- qr(P)
      P <- qr.Q(qrp)
      di <- ifelse( diag(qr.R(qrp))<0, -1, 1)
      P <- P * rep(di, each=nrow(P))
      # zapsmall( crossprod(P) ) # prove P'P orthogonal
    } else {
      eig <- sqrt(colSums(C^2))
    }
    
    if(verbose >  0){
      cat("EM Iter:", emiter, 
          " R2:", 1-(sum(x1*x1,na.rm=TRUE)/TotalSS), 
          "\n")
    }

    if( max(abs(P0-P)) < tol ) break # emiter
    
  } # Done finding PCs
  
  # output
  scores <- C
  loadings <- P
  
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
  
  # calculate R2 using non-missing cells
  R2cum <- rep(NA, ncomp)
  for(ii in 1:ncomp){
    xhat.ii <- tcrossprod(sweep( scores[,1:ii,drop=FALSE], 2, eig[1:ii], "*") ,
                          loadings[,1:ii,drop=FALSE])
    xhat.ii[x.miss] <- NA
    R2cum[ii] <- sum(xhat.ii * xhat.ii, na.rm=TRUE) / TotalSS
  }
  R2 <- c(R2cum[1], diff(R2cum))   # un-cumulate R2
  
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
              iter=emiter, 
              center=cmeans, scale=csds
  )
  return(out)
}


solve_weighted <- function(A, b, wt) {
  # Solve for x in wAx=wb with weight VECTOR wt
  # A matrix
  # b vector
  # wt vector of weights
  return(lm.wfit(A,b,wt)$coef)

  # same thing, but can have problems with singularities
  # wt <- sqrt(wt)
  # A <- wt * A # same as diag(wt) %*% A
  # b <- wt * b # same as diag(wt) %*% b
  # solve(crossprod(A), crossprod(A, b))
  # lm.fit(A,b)$coef

  # regress y on X, no intercept
  # coef(lm(y~ -1 + X))
  # coef(lm.fit(X,y)) is same as:
  # solve(t(X)%*%X) %*% t(X)%*%y
  
  #coef(lm(x[i,] ~ -1 + P))
  #lm.fit(P, x[i,])$coef
  #splom(cbind(x[i,],P), type=c("p","r")) # bottom row
}

