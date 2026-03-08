# Burr type XII distribution

Density, distribution function, quantile function and random generation
for the Burr type XII distribution on \\x\>0\\.

## Usage

``` r
dburr12(x, c, k, s, log = FALSE)

pburr12(q, c, k, s, lower.tail = TRUE, log.p = FALSE)

qburr12(p, c, k, s, lower.tail = TRUE, log.p = FALSE)

rburr12(n, c, k, s)
```

## Arguments

- x, q:

  Numeric vector of quantiles (\>0).

- c, k, s:

  Positive parameters.

- log:

  Logical; if TRUE return log-density.

- lower.tail:

  Logical; if TRUE (default) probabilities are \\P\[X \le x\]\\.

- log.p:

  Logical; if TRUE return log-probabilities.

- p:

  Numeric vector of probabilities between 0 and 1.

- n:

  Number of observations.

## Details

This parameterization requires `k > 1`.
