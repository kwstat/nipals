# test-uscrime.R

# see setup-expectations.R for expect_aligned

## ---------------------------------------------------------------------------

# This data is from the documentation for the SAS PCA procedure
# http://documentation.sas.com/?docsetId=casstat&docsetVersion=8.11&docsetTarget=viyastat_pca_gettingstarted.htm&locale=en
# Also HPPRINCOMP
# https://documentation.sas.com/doc/en/pgmsascdc/v_007/statug/statug_hpprincomp_gettingstarted.htm

data(uscrime)
# SAS discards rows with missinng values
dat <- uscrime[complete.cases(uscrime), ]
rownames(dat) <- dat[ ,1]
dat <- as.matrix(dat[ , -1])

sas_mean <- c(7.51667, 26.07500, 127.55625, 214.58750, 1316.37917, 2696.88542, 383.97917)

sas_sd <- c(3.93059, 10.81304, 88.49374, 100.64360, 423.31261, 714.75023, 194.37033)

sas_eig <- c(4.045824, 1.264030, 0.747500, 0.326325, 0.265207, 0.228364, 0.122750)

sas_loadings <- matrix(c(
  -0.30289,	0.61893,	0.17353,	-0.23308,	-0.54896,	0.26371,	-0.26428,
  -0.43410,	0.17053,	-0.23539,	0.06540,	-0.18075,	-0.78232,	0.27946,
  -0.39705,	-0.04713,	0.49208,	-0.57470,	0.50808,	-0.09452,	0.02497,
  -0.39622,	0.35142,	-0.05343,	0.61743,	0.51525,	0.17395,	-0.19921,
  -0.44164,	-0.20861,	-0.22454,	-0.02750,	-0.11273,	0.52340,	0.65085,
  -0.35634,	-0.40570,	-0.53681,	-0.23231,	-0.02172,	0.04085,	-0.60346,
  -0.28834,	-0.50400,	0.57524,	0.41853,	-0.35939,	-0.06024,	-0.15487),
  nrow=7, byrow=TRUE)
row.names(sas_loadings) <- c("murder", "rape", "robbery", "assault",
                             "burglary", "larceny", "autotheft")
colnames(sas_loadings) <- paste0("PC",1:7)

## ---------------------------------------------------------------------------

test_that("NIPALS decomposition of uscrime data matches SAS results", {
  
  m1 <- nipals(dat) # complete-data method
  
  m2 <- nipals(dat, force.na=TRUE) # use the missing-data method

  expect_equal(as.vector(m1$center), sas_mean)
  expect_equal(as.vector(m1$scale), sas_sd, tol=1e-5)
  
  # SAS eigenvalues reported are squared singular values divided by (n-1)
  expect_equal((m1$eig^2) / (nrow(dat)-1) , sas_eig, tolerance = 1e-4)
  expect_equal((m2$eig^2) / (nrow(dat)-1) , sas_eig, tolerance = 1e-4)
  
  expect_aligned(m1$loadings, sas_loadings, tol=1e-3)
  expect_aligned(m2$loadings, sas_loadings, tol=1e-3)

})

test_that("SVD, NIPALS and EMPCA results match for complete uscrime data", {

  m1 <- nipals::nipals(dat) # complete-data method
  m1s <- svd(scale(dat))

  # NIPALS vs SVD
  
  # eigenvalues
  expect_equal(m1$eig, m1s$d, tol = 1e-3)
  # scores
  expect_aligned(m1s$u, m1$scores)
  # loadings
  expect_aligned(m1s$v, m1$loadings)

  # EMPCA vs SVD
  
  m1e <- empca(dat)
  expect_equal(m1$eig, m1e$eig, tol=1e-3)
  expect_aligned(m1s$u, m1e$scores)  
  expect_aligned(m1s$v, m1e$loadings)
  
})

