# Derived Burr-like distribution

Density, distribution function, quantile function and random generation
for a derived Burr-like distribution with CDF \$\$F(x) = \left(1 +
(T/x)^a \exp((T-x)/b)\right)^{-1},\quad x\>0.\$\$

## Usage

``` r
dburrNEW(x, a, b, T, log = FALSE)

pburrNEW(q, a, b, T, lower.tail = TRUE, log.p = FALSE)

qburrNEW(
  p,
  a,
  b,
  T,
  lower.tail = TRUE,
  log.p = FALSE,
  tol = 1e-10,
  maxiter = 200
)

rburrNEW(n, a, b, T)
```

## Arguments

- x, q:

  Numeric vector of quantiles (\>0).

- a, b, T:

  Positive parameters.

- log:

  Logical; if TRUE return log-density.

- lower.tail:

  Logical; if TRUE (default) probabilities are \\P\[X \le x\]\\.

- log.p:

  Logical; if TRUE return log-probabilities.

- p:

  Numeric vector of probabilities between 0 and 1.

- tol:

  Numerical tolerance passed to
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html).

- maxiter:

  Maximum iterations passed to
  [`uniroot`](https://rdrr.io/r/stats/uniroot.html).

- n:

  Number of observations.
