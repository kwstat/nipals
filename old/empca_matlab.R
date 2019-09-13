empca_matlab <- function(x, w, ncomp=min(nrow(x), ncol(x)),
                  center=TRUE, scale=TRUE,
                  maxiter=20,
                  seed=1,
                  tol=1e-6,
                  fitted=FALSE,
                  verbose=FALSE) {

  x <- as.matrix(x) # in case it is a data.frame
  nvar <- nvar <- ncol(x)
  nobs <- nobs <- nrow(x)
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
  
  # position of NAs
  x.miss <- is.na(x)
  has.na <- any(x.miss)
  # Change missing values to 0 with 0 weight
  x[x.miss] <- 0
  w[x.miss] <- 0

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
  #P <- P[1:nvar, 1:ncomp]
  # colSums(P^2)           # each column is unit length
  # zapsmall( t(P) %*% P ) # P is orthonormal in columns

  # ported from matlab empca_w.m
  normvec <- function(x) x/sqrt(sum(x^2))
  P <- matrix(NA, nrow=nobs, ncol=ncomp)
  C <- matrix(NA, nrow=nvar, ncol=ncomp)
  for(h in 1:ncomp){
    # random dir
    P[,h] <- normvec( runif(nobs) )

    for(iter in 1:maxiter) {
      P0 <- P[,h]
      #browser()
      # E step
      C[,h] <- t(x) %*% P[,h, drop=FALSE]
      # M step
      # column h of C times each column of t(W)
      CW <- t(w) * C[,h]
      P[,h] <- rowSums(x * t(CW)) / t(CW) %*% C[,h,drop=FALSE]
      P[,h] <- normvec(P[,h])


      #browser()
      # Renormalize and re-orthogonalize columns of P via GramSchmidt
      # P <- pracma::gramSchmidt(P)$Q
      # need to use native qr() instead
      # https://stackoverflow.com/questions/15584221/gram-schmidt-with-r
      # P <- qr.Q(qr(P))
      # check P has orthonormal columns
      # zapsmall( t(P) %*% P )

      if( max(abs(P0-P[,h])) < tol ) break
    }
    # deflate X
    if(verbose)
      cat("h: ", h, "iter: ", iter, "\n")
    x <- x - P[,h] %*% t(C[,h]) # deflate matlab
  }

  # why should it be that P is orthonormal?
  # It's not when there are unequal weights
  # zapsmall(t(P)%*%P, 3)
  
  eig = sqrt(colSums(C^2))
  V = apply(C, 2, normvec)

  # P = 'scores' from nipals
  # sv = eigenvalues
  # v = 'loadings' from nipals
  
  # calculate xhat=fitted
  
  # output
  scores <- P
  rownames(scores) <- rownames(x)
  colnames(scores) <- paste("PC", 1:ncol(scores), sep="")
  loadings <- V
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

