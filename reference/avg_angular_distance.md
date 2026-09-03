# Average angular distance between two rotation matrices

For matrices A and B, calculate the angle between the column vectors of
A and the corresponding column vectors of B. Then average the angles.

## Usage

``` r
avg_angular_distance(A, B)
```

## Arguments

- A:

  Matrix

- B:

  Matrix

## Value

A single floating point number, in radians.

## Details

The results of the singular value decomposition X=USV' are unique, but
only up to a change of sign for columns of U, which indicates that the
axis is flipped.

## References

None

## Author

Kevin Wright

## Examples

``` r
# Example from https://math.stackexchange.com/questions/2113634/
rot1 <- matrix(c(-0.956395958, -0.292073218, 0.000084963,
                 0.292073230, -0.956395931, 0.000227268,
                 0.000014880, 0.000242173, 0.999999971),
               ncol=3, byrow=TRUE)
rot2 <- matrix(c(-0.956227882, -0.292623029, -0.000021887,
                 0.292623030, -0.956227882, -0.000024473,
                 -0.000013768, -0.000029806, 0.999999999),
               ncol=3, byrow=TRUE)
avg_angular_distance(rot1, rot2) # .0004950387
#> [1] 0.0004950387
```
