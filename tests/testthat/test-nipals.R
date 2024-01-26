# test-nipals.R

# see setup-expectations.R for expect_aligned

# Note, I looked at the python repository:
# https://github.com/fredrikw/python-nipals/blob/master/tests/test_nipals.py
# I did not understand the tests and so did not incorporate the yarn
# data or oliveoil data.

## ---------------------------------------------------------------------------

# Create data for tests

B_full <- matrix(c(50, 67, 90, 98, 120,
                   55, 71, 93, 102, 129,
                   65, 76, 95, 105, 134,
                   50, 80, 102, 130, 138,
                   60, 82, 97, 135, 151,
                   65, 89, 106, 137, 153,
                   75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
rownames(B_full) <- c("G1","G2","G3","G4","G5","G6","G7")
colnames(B_full) <- c("E1","E2","E3","E4","E5")

# Same data with 2 NAs
B_miss <- matrix(c(NA, 67, 90, 98, 120,
                   NA, 71, 93, 102, 129,
                   65, 76, 95, 105, 134,
                   50, 80, 102, 130, 138,
                   60, 82, 97, 135, 151,
                   65, 89, 106, 137, 153,
                   75, 95, 117, 133, 155), ncol=5, byrow=TRUE)
rownames(B_miss) <- c("G1","G2","G3","G4","G5","G6","G7")
colnames(B_miss) <- c("E1","E2","E3","E4","E5")

## ---------------------------------------------------------------------------

test_that("For complete-data matrix, fitted values match original data", {

  m1 <- nipals(B_full, center=TRUE, scale=TRUE, fitted=TRUE, tol=1e-8)
  expect_equal( round( m1$fitted, 3), B_full )

  m2 <- nipals(B_full, center=TRUE, scale=FALSE, tol=1e-8, fitted=TRUE)
  expect_equal( round( m2$fitted, 3), B_full )
  
  m3 <- nipals(B_full, center=FALSE, scale=FALSE, tol=1e-8, fitted=TRUE)
  expect_equal( round( m3$fitted, 3), B_full )
  
  m4 <- nipals(B_full, center=FALSE, scale=TRUE, tol=1e-8, fitted=TRUE)
  expect_equal( round( m4$fitted, 3), B_full )
    
})

test_that("Estimate fitted values for NAs", {
  # Published examples of NIPALS with missing values are almost nonexistant.
  # Here is an example with XLSTAT (details on methodology are limited)
  # 
  # Circa 2017 XLSTAT estimated the missing values below as
  # https://help.xlstat.com/customer/en/portal/articles/2062415-missing-data-imputation-using-nipals-in-excel?b_id=9283
  #   1365.236, 88.600, 175.798, 1051.698, 432.470, 168.554
  #
  # In 2019, XLSTAT estimates
  # https://help.xlstat.com/6497-missing-data-imputation-using-nipals-excel
  #   1398.829, 88.151, 178.313, 1060.534, 444.254, 169.686
  #
  auto <- data.frame(capacity = c(1396, 1721, 1580, 1769, 2068, 1769L),
                     power = c(90, 92, 83, 90, 88, 90L),
                     speed = c(174, 180, 170, 180, 180, 182L),
                     weight = c(850, 965, 970, 1080, 1135, 1060L),
                     length = c(369, 415, 395, 440, 446, 424L),
                     width = c(166, 169, 170, 169, 170, 168L))
  rownames(auto) <- c("Honda civic", "Renault 19", "Fiat Tipo",
                      "Peugeot 405", "Renault 21", "Citroen BX")
  # create missing values along the diagonal
  auto[1,1] <- auto[2,2] <- auto[3,3] <- auto[4,4] <- auto[5,5] <- auto[6,6] <- NA

  library(nipals)
  # These settings give results similar to XLstat 
  m1 <- nipals(auto, fitted=TRUE, gramschmidt=FALSE)
  expect_equal(diag(m1$fitted),
               c(1365.236, 88.600, 175.798, 1051.698, 432.470, 168.554),
               tol=1e-1)

  # Test Github issue #2
  expect_silent(nipals(auto[,-1], ncomp=1, fitted=TRUE))
})

test_that("Code coverage of nipals function arguments", {
  
  library(nipals)
  Bnarow = B_full
  Bnarow[1,] = NA
  expect_error(nipals(Bnarow))
  Bnacol = B_full
  Bnacol[,1] = NA
  expect_error(nipals(Bnacol))
  
  # ncomp
  m1 = nipals(B_miss, ncomp=1)
  
  # center/scale
  m1 = nipals(B_miss, center=FALSE, scale=FALSE)
  m1 = nipals(B_miss, center=TRUE, scale=FALSE)
  m1 = nipals(B_miss, center=TRUE, scale=TRUE)
  
  # maxiter
  expect_warning(nipals(B_full, maxiter=2))

  # tol
  m1 = nipals(B_miss, tol=1e-1, verbose=TRUE)
  m1 = nipals(B_miss, tol=1e-10, verbose=TRUE)

  # startcol
  m1 = nipals(B_full, startcol=0, verbose=TRUE)
  m1 = nipals(B_full, startcol=5, verbose=TRUE)

  # fitted
  expect_null(nipals(B_full)$fitted)
  expect_null(nipals(B_full, fitted=FALSE)$fitted)
  expect_false(is.null(nipals(B_full, fitted=TRUE)$fitted))

  # force.na
  m1 = nipals(B_full, force.na=FALSE)
  m1 = nipals(B_full, force.na=TRUE)
  
  # gramschmidt
  m1 = nipals(B_miss, gramschmidt=FALSE) # default
  round(crossprod(m1$loadings), 3)  # 1 on diagonal, but not identity
  m2 = nipals(B_miss, gramschmidt=TRUE)
  round(crossprod(m2$loadings), 3)  # should be identity 5x5
  
  # verbose
  m1 = nipals(B_full, verbose=TRUE)
  
})

test_that("Start column function", {
  # corn data from C.Majer
  corn <- structure(c(20.73, 23.58, 22.41, 19.97, 21.42, 24.48, 23.19, 25.73,
                      22.66, 23.93, 18.79, 18.56, 18.47, 18.33, 20.01, 19.03, 
                      19.36, 21.2, 19.17, 18.17, 20.41, 17.67, 17.89, 21.41,
                      18.81, 22.77, NA, 19.86, 19.43, NA, 17.28, 14.9, 16.52,
                      14.77, 19.35, 22.46, 24.61, NA, NA, 26.16, NA, 19.48,
                      NA, 20.72, 17.8, 18.55, 19.56, 20.01, 20.05, 17.18,
                      16.29, 17.41, 15.86, 15.7, 17.84, 24.1, 27.02, 26.76,
                      26, 26.02, 20.63, 20.37, 21.17, 21.55, 19.12),
                    .Dim = c(5L, 13L),
                    .Dimnames = list(c("G1", "G2", "G3", "G4", "G5"),
                                     c("E01", "E02", "E03", "E04", "E05",
                                       "E06", "E07", "E08", "E09", "E10",
                                       "E11", "E12", "E13")))
  # specify the startcol as a function
  expect_silent( nipals(corn, startcol=function(x) sum(abs(x), na.rm=TRUE)) )
  # using column with maximum variance fails
  # do NOT use this line...has problems when CRAN checks with --noLD
  # expect_warning( nipals(corn, startcol=function(x) var(x, na.rm=TRUE)) )
  
})

test_that("Calculate predictions from model", {
  
  # choose 75% of rows for training sample
  set.seed(42)
  ix <- sample(nrow(iris), nrow(iris)*0.75)
  iris.train <- iris[ix,1:4]
  iris.test <- iris[-ix,1:4]

  # Pre-computed predictions for reference
  p1ref <- structure(
    c(-2.1634, -2.45126, -2.53552, -2.32751, -2.41099, 
      -0.8063, -0.52704, -0.18264, 0.02092, -1.26167, 0.26757, -0.01858, 
      -0.31982, 0.10316, -0.10626, -0.09459, -0.02626, 0.03138, 0.02455, 
      0.0337), .Dim = 5:4, .Dimnames = list(c("2", "3", "7", "8", "9"),
                                            c("PC1", "PC2", "PC3", "PC4")))

  # Method 1: Assign a class to the nipals model
  m1 <- nipals(iris.train, startcol=2)
  class(m1) <- "princomp"
  p1 <- predict(m1, newdata=iris.test)

  # Method 2: Call stats:::predict.princomp
  # m2 <- nipals(iris.train[,1:4])
  # stats:::predict.princomp(m2, newdata=iris.test[,1:4])
    
  # # prcomp uses x$rotation, princomp uses x$loadings
  # m3 <- m1
  # m3$rotation <- m3$loadings
  # stats:::predict.prcomp(m3, newdata=iris.valid[,1:4])
  # 
  # # princomp 
  # m4 <- princomp(iris.train[,1:4])
  # predict(m4, newdata=iris.valid[,1:4])
  
})

test_that("avg_angular_distance", {
  #  https://math.stackexchange.com/questions/2113634/
  rot1 <- matrix(c(-0.956395958, -0.292073218, 0.000084963,
                   0.292073230, -0.956395931, 0.000227268,
                   0.000014880, 0.000242173, 0.999999971),
                 ncol=3, byrow=TRUE)
  rot2 <- matrix(c(-0.956227882, -0.292623029, -0.000021887,
                   0.292623030, -0.956227882, -0.000024473,
                   -0.000013768, -0.000029806, 0.999999999),
                 ncol=3, byrow=TRUE)
  expect_equal(avg_angular_distance(rot1, rot2), .0004950387)
})
