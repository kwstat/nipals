# empca.R

# Idea: Use a mixed model to calculate the biplot (similar to AMMI?).
# Perhaps this will 'weight' locations that have more data?
# Try subsetting a multi-rep dataset so that some locs have 1
# rep and some locs have more reps. Compare this to a weighted biplot of
# blues.

# NA values with weight 0
B0 <- matrix(c(50, 67, 90, 98, 120,
               55, 71, 93, 102, 129,
               65, 76, 95, 105, 134,
               50, 80, 102, 130, 138,
               60, 82, 97, 135, 151,
               65, 89, 106, 137, 153,
               75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
rownames(B0) <- c("G1","G2","G3","G4","G5","G6","G7")
colnames(B0) <- c("E1","E2","E3","E4","E5")
B0wt <- matrix(1, nrow=nrow(B0), ncol=ncol(B0))

# m1: NA values with weight 1
B1 <- B0
B1wt <- B0wt
B1[1,1] <- B1[2,1] <- NA
m1 <- empca(B1, B1wt)
# m2: NA values with weight 10
B2 <- B0
B2[1,1] = B2[2,1] = NA
B2wt <- B1wt
B2wt[1,1] <- B2wt[2,1] <- 10
m2 <- empca(B2, B2wt, verbose=TRUE)

# m3: Regular values with 0 weight
B3 <- B1
B3wt <- B1wt
B3wt[1,1] <- B3wt[2,1] <- 0
m3 <- empca(B3, B3wt) # oddly, this does not match B2
# m4: Outliers with 0 weight
B4 <- B1
B4wt <- B1wt
B4[1,1] <- B4[2,1] <- 200
B4wt[1,1] <- B4wt[2,1] <- 0
m4 <- empca(B4, B4wt) # oddly, this does not match B2


if(0){
  
  B1 <- matrix(c(50, 67, 90, 98, 120,
                55, 71, 93, 102, 129,
                65, 76, 95, 105, 134,
                50, 80, 102, 130, 138,
                60, 82, 97, 135, 151,
                65, 89, 106, 137, 153,
                75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
  rownames(B1) <- c("G1","G2","G3","G4","G5","G6","G7")
  colnames(B1) <- c("E1","E2","E3","E4","E5")
  B1wt <- matrix(1, nrow=nrow(B1), ncol=ncol(B1))

  B2 <- B1
  B2[1,1] = B2[2,1] = NA
  B2wt <- B1wt
  B2wt[1,1] <- B2wt[2,1] <- 0

  m1s <- svd(B1)
  m1s$u[,1] <- -1 * m1s$u[,1]
  m1s$v[,1] <- -1 * m1s$v[,1]
  
  m1n <- nipals::nipals(B1)
  #m1$loadings
  #m1n$loadings
  #m1$scores
  #m1n$scores

  biplot(m1s$u, m1s$v, main="SVD")
  biplot(m1n$scores, m1n$loadings, main="Nipals") # exactly the same as svd


  #m1p <- empca_python(x=B1, w=B1wt, ncomp=4)
  #m1mat <-  empca_matlab(x=B1, w=B1wt, ncomp=4, tol=1e-12, verbose=TRUE)
  m1e <- empca(x=B1, w=B1wt, ncomp=4, maxiter=100, fitted=TRUE)
  m1n$eig
  m1e$eig
  
  m1e$scores[,1:2] <- -1 * m1e$scores[,1:2]
  m1e$loadings[,1:2] <- -1 * m1e$loadings[,1:2]
  biplot(m1e$scores, m1e$loadings, main="EMPCA") # same as nipals!
  
  zapsmall( crossprod( m1$scores ))
  zapsmall( crossprod( m1$loadings ))
  zapsmall( abs(m1n$scores[,1:4]) - abs(m1e$scores[,1:4]) ) # very similar
  zapsmall( abs(m1n$loadings[,1:4]) - abs(m1e$loadings[,1:4]) ) # very similar
  zapsmall (abs(m1mat$loadings) - abs(m1$loadings) )
  

  # now with missing values
  m2n <- nipals::nipals(B2)
  m2e <- empca(x=B2, w=B2wt, ncomp=5, maxiter=100, fitted=TRUE)
  m2n$eig
  m2e$eig
  m2n$R2
  m2e$R2
  m2n$loadings
  m2e$loadings
  m2n$scores
  m2e$scores
  biplot(m2n$scores, m2n$loadings, main="B2 - NIPALS")
  m2e$scores[,2] <- -1 * m2e$scores[,2]
  m2e$loadings[,2] <- -1 * m2e$loadings[,2]
  biplot(m2e$scores, m2e$loadings, main="B2 - EMPCA") # same as nipals!

  # Compare to Matlab
  round( m1e$scores, 3)   # P matches
  round( m1e$loadings, 3) # v matches
  round( m2e$scores, 3)
  round( m2e$loadings, 3)

  set.seed(43)
  Bbig <- matrix(rnorm(100*100), nrow=100)
  Bbig2 <- Bbig
  Bbig2[1,1] <- NA
  system.time(m3n <- nipals::nipals(Bbig2, maxiter=1000)) # 3.3 sec
  system.time(m3e <- empca(x=Bbig2, maxiter=100)) # 10.4 sec
  h(m3n$eig)
  h(m3e$eig)
  # high agreement in PCs, horizontal stripes are probably swaps of two PCs
  levelplot(abs(m3n$scores) - abs(m3e$scores), col.regions=gge::RedGrayBlue)
  levelplot(abs(m3n$loadings) - abs(m3e$loadings), col.regions=gge::RedGrayBlue)
}


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
    P0 <- P # previous P to check convergence

    # E step: find C[i,k] such that X[i,] = Sum_k: C[i,k] P[,k]
    for(i in 1:nobs) {
      # C[i,] <- t( solve_weighted(P, x[i,], w[i,]) )
      C[i,] <- t(lm.wfit(P, x[i,], w[i,])$coef)
    }

    # replace NA by 0
    C[is.na(C)] <- 0

    # M step: Calculate eigvectors in P
    x1 = x # reset for each EM iteration
    for( h in 1:ncomp ) {
      # Bailey eqn 21.
      # P[j,h] = sum_j (x[,j] * w[,j] * C[,h]) / sum_j (w[,j]*C[,h]*C[,h])
      # P[,h] = colSums(x * w *C[,h]) / colSums(w * C[,h] * C[,h])
      wC = w*C[,h] # column C[,h] times each column of w
      P[,h] <- colSums(x1 * wC) / colSums(wC * C[,h])
      # If a column of C is entirely NA, then P[,h] is NA. Should we just
      # set it 0 instead?
      if(all(is.na(P[,h]))) P[,h] <- 0
      # not necessary to unitize columns P[,h]
      # x1 <- x1 - P[,h] %*% t(C[,h]) # deflate, Bailey eqn 22.
      x1 <- x1 - outer(C[,h],P[,h])
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

      # perform gram-schmidt rotation on C
      qrc <- qr(C)
      C <- qr.Q(qrc) # C is now orthonormal with unit-length columns
      # flip sign of columns where diagonal elements of R are negative
      di <- diag(qr.R(qrc))
      di <- di/abs(di) # vector of +1 and -1
      C <- C * rep(di, each=nrow(C)) # flip sign of columns as needed

      # do we really need to do this, or can we get information from the
      # gram-schmidt rotation of C ???
      
      # gram-schmidt on P
      qrp <- qr(P)
      P <- qr.Q(qrp)
      di <- diag(qr.R(qrp))
      di <- di/abs(di)
      P <- P * rep(di, each=nrow(P))
      # zapsmall( t(P) %*% P ) # prove orthogonal
    } else {
      eig <- sqrt(colSums(C^2))
    }
    
    if(verbose){
      cat("EM Iter:", emiter, " R2:", 1-(sum(x1*x1,na.rm=TRUE)/TotalSS), "\n")
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
}

