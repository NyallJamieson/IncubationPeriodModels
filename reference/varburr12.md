# Variance of the Burr type XII distribution

Variance for the Burr type XII distribution (this parameterization
requires `k > 1`). The variance exists only if `k > 1 + 2/c`.

## Usage

``` r
varburr12(c, k, s)
```

## Arguments

- c, k, s:

  Positive parameters.

## Value

A numeric value (may be `Inf` if the variance does not exist).
