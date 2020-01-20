# empca_python.R

empca_python <- function(x, w, ncomp=min(nrow(x), ncol(x)),
                  center=TRUE, scale=TRUE,
                  maxiter=20,
                  seed=1,
                  tol=1e-6,
                  fitted=FALSE,
                  verbose=FALSE) {

  x <- as.matrix(x) # in case it is a data.frame
  nvar <- ncol(x)
  nobs <- nrow(x)
  x.orig <- x # Save x for replacing missing values

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
  
  # initialize outputs
  C <- matrix(NA, nrow=nobs, ncol=ncomp)
  R2cum <- rep(NA, length=ncomp)

  # position of NAs
  x.miss <- is.na(x)
  has.na <- any(x.miss)
  x[x.miss] <- 0

  # initial E step, P matrix with random orthonormal columns, nvar*ncomp
  # derived from ideas in pracma::randortho
  set.seed(seed)
  dd <- max(ncomp,nvar)
  ranmat <- matrix(rnorm(dd*dd), nrow=dd)
  qrmod <- qr(ranmat)
  di <- diag(qr.R(qrmod))
  # di <- di/abs(di) # di is vector of +1 and -1. Why needed???
  qmat <- qr.Q(qrmod)
  P <- qmat %*% diag( di/abs(di) ) %*% qmat
  # fixme - use identity matrix for testing purposes
  P <- diag(max(nvar, ncomp))
  P <- P[1:nvar, 1:ncomp]
  # colSums(P^2)           # each column is unit length
  # zapsmall( t(P) %*% P ) # P is orthonormal in columns

  for(iter in 1:maxiter){
    # E step: find C[i,k] such that X[i,] = Sum_k: C[i,k] P[,k]
    for(i in 1:nobs) {
      cat(i)
      C[i,] <- t( solve_weighted(P, x[i,], w[i,]) )
    }
    # replace NA by 0
    C[is.na(C)] <- 0
    # M step. Find P
    P <- calc_eigvec(x, w, C, ncomp=ncomp)
    
    # Calc R2, Rchi2
    xhat <- C %*% t(P)
    ssdev <- sum((xhat-x)^2, na.rm=TRUE)
    # check convergence
    if(verbose) cat("iter ", iter, "ssdev", ssdev, "R2", "Rchi2", "\n")   
  }

  ## # ported from matlab empca_w.m
  normvec <- function(x) x/sqrt(sum(x^2))

  eig = sqrt(colSums(C^2))
  V = apply(C, 2, normvec)
  
  # P = 'scores' from nipals
  # sv = eigenvalues
  # v = 'loadings' from nipals
  

  # output
  scores <- C
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste("PC", 1:ncol(scores), sep="")
  loadings <- P
  rownames(loadings) <- colnames(x)
  colnames(loadings) <- paste("PC", 1:ncol(loadings), sep="")
  
  out <- list(eig=eig,
              scores=scores,
              loadings=loadings,
              #fitted=xhat,
              ncomp=ncomp,
              #R2=R2,
              iter=iter, 
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

  # this can have problems with singularities
  wt <- sqrt(wt)
  A <- wt * A # same as diag(wt) %*% A
  b <- wt * b # same as diag(wt) %*% b
  #solve(crossprod(A), crossprod(A, b))
  lm.fit(A,b)$coef
}


calc_eigvec <- function(x, w, C, ncomp) {
  # Solve for eigvec P[i,] such that X[i,] = Sum_k: C[i,k] P[,k]

  # Calculate the eigenvectors one by one, up to ncomp
  
  # fixme: can we remove one loop via apply?
  P <- matrix(NA, nrow=ncol(x), ncol=ncomp)
  for(h in 1:ncomp) {
    ch <- C[, h]
    for(j in 1:ncol(x)) {
      wj <- w[, j]
      xj <- x[, j]
      cw <- ch * wj
      # P[j,h] = sum_j (x[,j] * w[,j] * C[,h]) / sum_j (C[,h]*C[,h]*w[,j])
      P[j,h]  <- sum(xj * cw) / sum(ch * cw)
    }
    
    x  <- x - outer(C[,h], P[,h])
  }

  # Renormalize and re-orthogonalize columns of P via GramSchmidt
  # P <- pracma::gramSchmidt(P)$Q
  # need to use native qr() instead
  # https://stackoverflow.com/questions/15584221/gram-schmidt-with-r
  P <- qr.Q(qr(P))
  # check P has orthonormal columns
  # zapsmall( t(P) %*% P )

  return(P)
}


if(FALSE) {
  # notes for solve_weighted
  # https://stats.stackexchange.com/questions/218146/weighted-least-square-weights-definition-r-lm-function-vs-mathbf-w-mathbf-a
  A1 <- as.matrix(mtcars[ , c("wt","hp","disp")])
  b1 <- as.matrix(mtcars$mpg)
  wt1 <- round(runif(nrow(A1)),2)
  # Usual way to do weighted regression
  lm(mpg ~ wt+hp+disp-1, data=mtcars, weights=wt1)
  # Hand-calculation: Pre-multiply by sqrt(diag(wt)), then solve()
  A2 <- sqrt(wt1) * A1
  b2 <- sqrt(wt1) * b1
  # because, t(A2) %*% A2 = crossprod(A2, A2) = crossprod(A2)
  # following 3 rows are the same
  # solve(t(A2) %*% A2, t(A2) %*% b2)
  # solve(crossprod(A2, A2), crossprod(A2, b2))
  solve(crossprod(A2), crossprod(A2, b2))
  lm.wfit(A1, b1, wt1)$coef
  lm.fit(A2, b2)$coef
  .lm.fit(A2, b2)$coef # do NOT use...moves singular terms to end
  solve_weighted(A1, b1, wt1)
}
