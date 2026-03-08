# Mean of the Burr type III distribution

Computes the mean of the Burr type III distribution on \\x\>0\\. The
mean is computed numerically from the survival function: \$\$E\[X\] =
\int_0^\infty (1 - F(x))\\dx.\$\$

## Usage

``` r
meanburr3(c, k, s, rel.tol = 1e-10)
```

## Arguments

- c, k, s:

  Positive parameters.

- rel.tol:

  Relative tolerance passed to
  [`integrate`](https://rdrr.io/r/stats/integrate.html).

## Value

A numeric value (may be `Inf` or `NA` if not finite / not computable).
