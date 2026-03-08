# Variance of the Burr type III distribution

Computes the variance of the Burr type III distribution on \\x\>0\\. The
variance is computed numerically from moments of the density.

## Usage

``` r
varburr3(c, k, s, rel.tol = 1e-10)
```

## Arguments

- c, k, s:

  Positive parameters.

- rel.tol:

  Relative tolerance passed to
  [`integrate`](https://rdrr.io/r/stats/integrate.html).

## Value

A numeric value (may be `Inf` or `NA` if not finite / not computable).
