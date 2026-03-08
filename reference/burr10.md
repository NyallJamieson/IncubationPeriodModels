# Burr type X distribution

Density, distribution function, quantile function and random generation
for the Burr type X distribution on \\x\>0\\.

## Usage

``` r
dburr10(x, c, k, log = FALSE)

pburr10(q, c, k, lower.tail = TRUE, log.p = FALSE)

qburr10(p, c, k, lower.tail = TRUE, log.p = FALSE)

rburr10(n, c, k)
```

## Arguments

- x, q:

  Numeric vector of quantiles (\>0).

- c:

  Positive parameter.

- k:

  Positive parameter.

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
