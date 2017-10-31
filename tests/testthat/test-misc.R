context("misc")

test_that("Code coverage of function arguments", {
  B <- matrix(c(50, 67, 90, 98, 120,
                55, 71, 93, 102, 129,
                65, 76, 95, 105, 134,
                50, 80, 102, 130, 138,
                60, 82, 97, 135, 151,
                65, 89, 106, 137, 153,
                75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
  rownames(B) <- c("G1","G2","G3","G4","G5","G6","G7")
  colnames(B) <- c("E1","E2","E3","E4","E5")
  
  library(nipals)
  Bnarow = B
  Bnarow[1,] = NA
  expect_error(nipals(Bnarow))
  Bnacol = B
  Bnacol[,1] = NA
  expect_error(nipals(Bnacol))
  
  B2 = B
  B2[1,1] = B2[2,1] = NA
  
  # ncomp
  m1 = nipals(B2, ncomp=1)
  
  # center/scale
  m1 = nipals(B2, center=FALSE, scale=FALSE)
  m1 = nipals(B2, center=TRUE, scale=FALSE)
  m1 = nipals(B2, center=TRUE, scale=TRUE)
  
  # maxiter
  expect_warning(nipals(B, maxiter=2))

  # tol
  m1 = nipals(B2, tol=1e-1, verbose=TRUE)
  m1 = nipals(B2, tol=1e-10, verbose=TRUE)

  # startcol
  m1 = nipals(B, startcol=0, verbose=TRUE)
  m1 = nipals(B, startcol=5, verbose=TRUE)

  # fitted
  expect_null(nipals(B)$fitted)
  expect_null(nipals(B, fitted=FALSE)$fitted)
  expect_false(is.null(nipals(B, fitted=TRUE)$fitted))

  # force.na
  m1 = nipals(B, force.na=FALSE)
  m1 = nipals(B, force.na=TRUE)
  
  # gramschmidt
  m1 = nipals(B2, gramschmidt=FALSE) # default
  round(crossprod(m1$loadings), 3)  # 1 on diagonal, but not identity
  m2 = nipals(B2, gramschmidt=TRUE)
  round(crossprod(m2$loadings), 3)  # should be identity 5x5
  
  # verbose
  m1 = nipals(B, verbose=TRUE)
  
})
