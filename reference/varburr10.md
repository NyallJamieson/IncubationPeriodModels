# Variance of the Burr type X distribution

Computes the variance of the Burr type X distribution on \\x\>0\\. The
variance is computed numerically from moments of the density.

## Usage

``` r
varburr10(c, k, rel.tol = 1e-10)
```

## Arguments

- c:

  Positive parameter.

- k:

  Positive parameter.

- rel.tol:

  Relative tolerance passed to
  [`integrate`](https://rdrr.io/r/stats/integrate.html).

## Value

A numeric value (may be `Inf` or `NA` if not finite / not computable).
