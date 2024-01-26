# setup-expectations.R

expect_aligned <- function(A, B, tol=.01){
  expect_lt( avg_angular_distance(A, B), tol)
}
