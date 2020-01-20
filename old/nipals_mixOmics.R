
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

m3 <- nipals(scale(B2), ncomp=5)

# note, mixOmics does not divide by (nr-1), but does take sqrt
m3$eig
# 4.8762447 2.0442446 1.0728240 0.2370526 0.1432616

## R> m3$eig^2
## [1] 23.77776272  4.17893615  1.15095130  0.05619393  0.02052388

sum(m3$eig^2) # 29.18437 # note: sum(scale(B2)^2, na.rm=TRUE) = 28

# P loadings, match ade4, plsdepot
m3$p

# T scores
m3$t # has eigen values swept out
m3$t* rep(m3$eig,each=7) # match

# svd singular values
# R> svd(scale(B))$d
# [1] 5.02051843 1.87932352 1.10817656 0.17225189 0.06936698

# ----- timing example -----
# 100 x 100
set.seed(43)
Bbig <- matrix(rnorm(100*100), nrow=100)
Bbig2 <- Bbig
Bbig2[1,1] <- NA
system.time(nipals(scale(Bbig2), ncomp=1)) # Only 1 factor !
#   user  system elapsed 
#   0.06    0.00    0.06 
system.time(res <- nipals(scale(Bbig2), ncomp=100)) # 100 factors
#   user  system elapsed 
#  19.91    0.50   26.27 
system.time(nipals(scale(Bbig), ncomp=100)) # 100 factors

nipals = function (X, ncomp = 1, reconst = FALSE, max.iter = 500, tol = 1e-09) {
    
    #-- X matrix
    if (is.data.frame(X))
    X = as.matrix(X)
    
    if (!is.matrix(X) || is.character(X))
    stop("'X' must be a numeric matrix.", call. = FALSE)
    
    if (any(apply(X, 1, is.infinite)))
    stop("infinite values in 'X'.", call. = FALSE)
    
    nc = ncol(X)
    nr = nrow(X)
    #-- put a names on the rows and columns of X --#
    X.names = colnames(X)
    if (is.null(X.names))
    X.names = paste("V", 1:ncol(X), sep = "")
    
    ind.names = rownames(X)
    if (is.null(ind.names))
    ind.names = 1:nrow(X)
    
    #-- ncomp
    if (is.null(ncomp) || !is.numeric(ncomp) || ncomp < 1 || !is.finite(ncomp))
    stop("invalid value for 'ncomp'.", call. = FALSE)
    
    #-- reconst
    if (!is.logical(reconst))
    stop("'reconst' must be a logical constant (TRUE or FALSE).",
    call. = FALSE)
    
    #-- max.iter
    if (is.null(max.iter) || max.iter < 1 || !is.finite(max.iter))
    stop("invalid value for 'max.iter'.", call. = FALSE)
    
    max.iter = round(max.iter)
    
    #-- tol
    if (is.null(tol) || tol < 0 || !is.finite(tol))
    stop("invalid value for 'tol'.", call. = FALSE)
    
    #-- pca approach -----------------------------------------------------------#
    #---------------------------------------------------------------------------#
    
    
    #-- initialisation des matrices --#
    p = matrix(nrow = nc, ncol = ncomp)
    t.mat = matrix(nrow = nr, ncol = ncomp)
    eig = vector("numeric", length = ncomp)
    nc.ones = rep(1, nc)
    nr.ones = rep(1, nr)
    is.na.X = is.na(X)
    na.X = FALSE
    if (any(is.na.X)) na.X = TRUE
    
    #-- boucle sur h --#
    for (h in 1:ncomp) {

      # kw change start column to 1
      #th = X[, which.max(apply(X, 2, var, na.rm = TRUE))]
      th = X[, 1]
      if (any(is.na(th))) th[is.na(th)] = 0
      ph.old = rep(1 / sqrt(nc), nc)
      ph.new = vector("numeric", length = nc)
      iter = 1
      diff = 1
      
      if (na.X) {
        X.aux = X
        X.aux[is.na.X] = 0
      }
        
      while (diff > tol & iter <= max.iter) {
        if (na.X) {
          ph.new = crossprod(X.aux, th)
          Th = drop(th) %o% nc.ones
          Th[is.na.X] = 0
          th.cross = crossprod(Th)
          ph.new = ph.new / diag(th.cross)
        } else {
          ph.new = crossprod(X, th) / drop(crossprod(th))
        }
        
        ph.new = ph.new / drop(sqrt(crossprod(ph.new)))
        
        if (na.X) {
          th = X.aux %*% ph.new
          P = drop(ph.new) %o% nr.ones
          P[t(is.na.X)] = 0
          ph.cross = crossprod(P)
          th = th / diag(ph.cross)  # <--------------
        } else {
          th = X %*% ph.new / drop(crossprod(ph.new))
        }
        
        diff = drop(sum((ph.new - ph.old)^2, na.rm = TRUE))
        ph.old = ph.new
        iter = iter + 1
      }
      
      if (iter > max.iter)
        warning(paste("Maximum number of iterations reached for comp.", h))
      
      X = X - th %*% t(ph.new)
      p[, h] = ph.new
      t.mat[, h] = th
      eig[h] = sum(th * th, na.rm = TRUE)
    }
  
  eig = sqrt(eig)
  t.mat = scale(t.mat, center = FALSE, scale = eig)
  attr(t.mat, "scaled:scale") = NULL
  result = list(eig = eig, p = p, t = t.mat)
  
  if (reconst) {
    X.hat = matrix(0, nrow = nr, ncol = nc)
    
    for (h in 1:ncomp) {
      X.hat = X.hat + eig[h] * t.mat[, h] %*% t(p[, h])
    }
        
    colnames(X.hat) = colnames(X)
    rownames(X.hat) = rownames(X)
    result$rec = X.hat
  }
  
  return(invisible(result))
}

