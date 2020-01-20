<<<<<<< HEAD
=======
# setup-expectations.R
>>>>>>> empca

expect_aligned <- function(A, B, tol=.01){

  # The results of SVD are unique, but only up to a change of sign
  # for columns of U, which indicates that the axis is flipped.
  # We need a way to see if two rotation matrices are almost aligned,
<<<<<<< HEAD
  # allowed for axis flips.
=======
  # allowing for axis flips.
>>>>>>> empca
  
  # Calculate the angle between each pair of columns from A and B,
  # and then average.
  # https://math.stackexchange.com/questions/2113634/

  dotprod <- colSums(A*B) / sqrt( colSums(A*A) * colSums(B*B) )
  # use abs() so that vectors which point in nearly opposite
  # directions ( dotprod= -1 ) can be considered as if pointing
  # in the same direction.
  # dotprod could be slightly larger than 1.0, so subtract eps.
  theta <- acos( abs(dotprod) - .Machine$double.eps )
  # average radian angle between pair-wise columns of A and B
  avg_angle <- mean(theta)
  expect_lt(avg_angle, tol)
}

