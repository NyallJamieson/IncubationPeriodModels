# Mean of the derived Burr-like distribution

Computes the mean numerically from the survival function: \$\$E\[X\] =
\int_0^\infty (1 - F(x))\\dx.\$\$

## Usage

``` r
meanburrNEW(a, b, T, rel.tol = 1e-10)
```

## Arguments

- a, b, T:

  Positive parameters.

- rel.tol:

  Relative tolerance passed to
  [`integrate`](https://rdrr.io/r/stats/integrate.html).

## Value

A numeric value (may be `Inf` or `NA` if not finite / not computable).
