# Variance of the derived Burr-like distribution

Computes the variance numerically from moments of the density.

## Usage

``` r
varburrNEW(a, b, T, rel.tol = 1e-10)
```

## Arguments

- a, b, T:

  Positive parameters.

- rel.tol:

  Relative tolerance passed to
  [`integrate`](https://rdrr.io/r/stats/integrate.html).

## Value

A numeric value (may be `Inf` or `NA` if not finite / not computable).
