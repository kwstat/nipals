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

m1 <- nipals(B2,nf=5)

m1$eig
sum(m1$eig^2) # 16.22

# P loadings
m1$c1 # P loadings
round(t(m1$c1) %*% m1$c1, 3) # not orthogonal

# T scores
m1$li # T scores

# ----- timing example -----
set.seed(43)
Bbig <- matrix(rnorm(100*100), nrow=100)
Bbig[1,1] <- NA
system.time(nipals(Bbig, nf=1)) # SLOW. Only 1 factor!
##  user  system elapsed 
## 39.98    0.00   40.03 







"nipals" <- function(df, nf=2, rec=FALSE,niter=100, tol = 1e-9){
  # df est un data frame contenant eventuellement des valeurs manquantes (NA)
  # nf nombre de facteurs a conserver
  # rec, si rec=T, la reconstitution des donnees sur les nf premiers axes est realisee
  # **********************************************************************************
  # df is a data frame which can contain missing values (NA)
  # nf number of axes to keep
  # rec, if rec=T, data recontsitution is performed with the nf first axes
  # n.max.iter= maximum number of iterations
  
  df <- data.frame(df)
  tol <- 1e-9 # tol pour la convergence
  nc <- ncol(df)
  nr <- nrow(df)
  nr.na <- apply(df, 2, function(x) sum(!is.na(x)))
  if (rec)
    x <- list(li=matrix(0,nr,nf),c1=matrix(0,nc,nf),co=matrix(0,nc,nf),
            eig=rep(0,nf),nb=rep(0,nf),rec=matrix(0,nr,nc))
  else
    x <- list(li=matrix(0,nr,nf),c1=matrix(0,nc,nf),co=matrix(0,nc,nf),
            eig=rep(0,nf),nb=rep(0,nf))
  row.names(x$c1) <- names(df)
  row.names(x$co) <- names(df)
  row.names(x$li) <- row.names(df)
  
  #X <- scale(df, center=T, scale=T, na.rm=TRUE)
  cmeans <- colMeans(df, na.rm=TRUE)
  # kw: take off the shrinkage factor
  csd <- apply(df, 2, sd, na.rm=TRUE) #* (nr.na - 1) / nr.na
  X <- sweep(sweep(df, 2, cmeans, "-"), 2, csd, "/")
  x$tab <- X
  for (h in 1:nf) {
    th <- X[,1]
    ph1 <- rep(1/sqrt(nc),nc)
    ph2 <- rep(1/sqrt(nc),nc)
    diff <- rep(1,nc)
    nb <- 0
    while (sum(diff^2, na.rm=TRUE)>tol & nb<=niter) {
      for (i in 1:nc) {
        the <- th[!is.na(X[,i])]
        ph2[i] <- sum(X[,i]*th, na.rm=TRUE)/sum(the*the,na.rm=TRUE)
      }
      ph2 <- ph2/sqrt(sum(ph2*ph2,na.rm=TRUE))
      for (i in 1:nr) {
        ph2e <- ph2[!is.na(X[i,])]
        th[i] <- sum(X[i,]*ph2, na.rm=TRUE)/sum(ph2e*ph2e,na.rm=TRUE)
      }
      diff <- ph2-ph1
      ph1 <- ph2
      nb <- nb+1
    }
    if(nb>niter) stop(paste("Maximum number of iterations reached for axis", h))
    X <- X-th%*%t(ph1)
    x$nb[h] <- nb # nombre d'iterations (number of iterations)
    x$li[,h] <- th # coordonnees des lignes (row coordinates)
    x$c1[,h] <- ph1 # coordonnees des colonnes de variance unit' (columns coordinates of unit variance)
    x$eig[h] <- sum(th*th,na.rm=TRUE)/(nr-1) # valeurs propres (pseudo-eigenvalues)
    x$co[,h] <- x$c1[,h]*sqrt(x$eig[h]) # coord. col. de variance lambda (column coordinates of variance lambda)
    
  }
  if (rec) {
    for (h in 1:nf) {
      x$rec <- x$rec+x$li[,h]%*%t(x$c1[,h]) # tableau reconstitue (reconstitued data)
    }
  }
  if (rec){
    x$rec=as.data.frame(x$rec)
    names(x$rec) <- names (df)
    row.names(x$rec) <- row.names(df)
  }
  
  x$call <- match.call()
  x$nf <- nf
  class(x) <- "nipals"
  if(any(diff(x$eig)>0)) warning("Eigenvalues are not in decreasing order. Results of the analysis could be problematics")
  return(x)
}
