# Mean of the Burr type XII distribution

Mean for the Burr type XII distribution (this parameterization requires
`k > 1`). The mean exists only if `k > 1 + 1/c`.

## Usage

``` r
meanburr12(c, k, s)
```

## Arguments

- c, k, s:

  Positive parameters.

## Value

A numeric value (may be `Inf` if the mean does not exist).
