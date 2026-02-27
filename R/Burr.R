#' Burr distributions (Type III, X, XII) and a derived Burr-like distribution
#'
#' This package implements density, distribution, quantile, and random generation
#' functions for several Burr distributions and a derived Burr-like distribution.
#'
#' Included families:
#' \itemize{
#'   \item Burr Type III: \code{dburr3}, \code{pburr3}, \code{qburr3}, \code{rburr3}
#'   \item Burr Type X: \code{dburr10}, \code{pburr10}, \code{qburr10}, \code{rburr10}
#'   \item Burr Type XII: \code{dburr12}, \code{pburr12}, \code{qburr12}, \code{rburr12}
#'   \item Derived distribution: \code{dburrNEW}, \code{pburrNEW}, \code{qburrNEW}, \code{rburrNEW}
#' }
#'
#' All distributions are supported on \eqn{x > 0}.
#'
#' @name burrdist
#' @keywords distribution
NULL

# -----------------------------
# Helpers (internal)
# -----------------------------
.check_pos <- function(...) {
  pars <- list(...)
  for (nm in names(pars)) {
    v <- pars[[nm]]
    if (any(!is.finite(v))) stop(nm, " must be finite.", call. = FALSE)
    if (any(v <= 0)) stop(nm, " must be > 0.", call. = FALSE)
  }
  invisible(TRUE)
}

.as_prob <- function(p, log.p) {
  if (log.p) exp(p) else p
}

# -----------------------------
# Burr Type III topic
# -----------------------------

#' Burr type III distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Burr type III distribution on \eqn{x>0}.
#'
#' @param x,q Numeric vector of quantiles (>0).
#' @param p Numeric vector of probabilities between 0 and 1.
#' @param n Number of observations.
#' @param c,k,s Positive parameters.
#' @param log Logical; if TRUE return log-density.
#' @param lower.tail Logical; if TRUE (default) probabilities are \eqn{P[X \le x]}.
#' @param log.p Logical; if TRUE return log-probabilities.
#' @return A numeric vector.
#' @name burr3
NULL

#' @rdname burr3
#' @export
dburr3 <- function(x, c, k, s, log = FALSE) {
  .check_pos(c = c, k = k, s = s)
  x <- as.numeric(x)

  out <- rep.int(if (log) -Inf else 0, length(x))
  ok <- is.finite(x) & x > 0
  if (!any(ok)) return(out)

  z <- x[ok] / s
  logdens <- log(c) + log(k) - log(s) + (-c - 1) * log(z) + (-k - 1) * log1p(z^(-c))
  out[ok] <- if (log) logdens else exp(logdens)
  out
}

#' @rdname burr3
#' @export
pburr3 <- function(q, c, k, s, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(c = c, k = k, s = s)
  q <- as.numeric(q)

  p <- numeric(length(q))
  p[!is.finite(q)] <- NA_real_
  p[q <= 0] <- 0

  ok <- is.finite(q) & q > 0
  if (any(ok)) {
    z <- q[ok] / s
    p[ok] <- exp(-k * log1p(z^(-c)))
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname burr3
#' @export
qburr3 <- function(p, c, k, s, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(c = c, k = k, s = s)
  p <- .as_prob(p, log.p)
  if (!lower.tail) p <- 1 - p

  out <- rep.int(NA_real_, length(p))
  bad <- !is.finite(p) | p < 0 | p > 1
  out[bad] <- NaN
  out[p == 0] <- 0
  out[p == 1] <- Inf

  ok <- is.finite(p) & p > 0 & p < 1
  if (any(ok)) {
    out[ok] <- s * (p[ok]^(-1 / k) - 1)^(-1 / c)
  }
  out
}

#' @rdname burr3
#' @export
rburr3 <- function(n, c, k, s) {
  .check_pos(c = c, k = k, s = s)
  n <- as.integer(n)
  if (length(n) != 1 || n < 0) stop("n must be a nonnegative integer.", call. = FALSE)
  qburr3(stats::runif(n), c = c, k = k, s = s)
}

# -----------------------------
# Burr Type X topic
# -----------------------------

#' Burr type X distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Burr type X distribution on \eqn{x>0}.
#'
#' @inheritParams burr3
#' @param c Positive parameter.
#' @param k Positive parameter.
#' @name burr10
NULL

#' @rdname burr10
#' @export
dburr10 <- function(x, c, k, log = FALSE) {
  .check_pos(c = c, k = k)
  x <- as.numeric(x)

  out <- rep.int(if (log) -Inf else 0, length(x))
  ok <- is.finite(x) & x > 0
  if (!any(ok)) return(out)

  z <- (x[ok] / c)^2
  logdens <- log(2) + log(k) + log(x[ok]) - 2 * log(c) - z + (k - 1) * log1p(-exp(-z))
  out[ok] <- if (log) logdens else exp(logdens)
  out
}

#' @rdname burr10
#' @export
pburr10 <- function(q, c, k, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(c = c, k = k)
  q <- as.numeric(q)

  p <- numeric(length(q))
  p[!is.finite(q)] <- NA_real_
  p[q <= 0] <- 0

  ok <- is.finite(q) & q > 0
  if (any(ok)) {
    z <- (q[ok] / c)^2
    p[ok] <- (1 - exp(-z))^k
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname burr10
#' @export
qburr10 <- function(p, c, k, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(c = c, k = k)
  p <- .as_prob(p, log.p)
  if (!lower.tail) p <- 1 - p

  out <- rep.int(NA_real_, length(p))
  bad <- !is.finite(p) | p < 0 | p > 1
  out[bad] <- NaN
  out[p == 0] <- 0
  out[p == 1] <- Inf

  ok <- is.finite(p) & p > 0 & p < 1
  if (any(ok)) {
    out[ok] <- c * sqrt(-log1p(-p[ok]^(1 / k)))
  }
  out
}

#' @rdname burr10
#' @export
rburr10 <- function(n, c, k) {
  .check_pos(c = c, k = k)
  n <- as.integer(n)
  if (length(n) != 1 || n < 0) stop("n must be a nonnegative integer.", call. = FALSE)
  qburr10(stats::runif(n), c = c, k = k)
}

# -----------------------------
# Burr Type XII topic
# -----------------------------

#' Burr type XII distribution
#'
#' Density, distribution function, quantile function and random generation
#' for the Burr type XII distribution on \eqn{x>0}.
#'
#' This parameterization requires \code{k > 1}.
#'
#' @inheritParams burr3
#' @name burr12
NULL

#' @rdname burr12
#' @export
dburr12 <- function(x, c, k, s, log = FALSE) {
  .check_pos(c = c, s = s)
  if (any(k <= 1)) stop("k must be > 1 for this parameterization.", call. = FALSE)

  x <- as.numeric(x)
  out <- rep.int(if (log) -Inf else 0, length(x))
  ok <- is.finite(x) & x > 0
  if (!any(ok)) return(out)

  z <- x[ok] / s
  logdens <- log(c) + log(k - 1) - log(s) + (c - 1) * log(z) - k * log1p(z^c)
  out[ok] <- if (log) logdens else exp(logdens)
  out
}

#' @rdname burr12
#' @export
pburr12 <- function(q, c, k, s, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(c = c, s = s)
  if (any(k <= 1)) stop("k must be > 1 for this parameterization.", call. = FALSE)

  q <- as.numeric(q)
  p <- numeric(length(q))
  p[!is.finite(q)] <- NA_real_
  p[q <= 0] <- 0

  ok <- is.finite(q) & q > 0
  if (any(ok)) {
    z <- (q[ok] / s)^c
    p[ok] <- 1 - exp((1 - k) * log1p(z))
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname burr12
#' @export
qburr12 <- function(p, c, k, s, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(c = c, s = s)
  if (any(k <= 1)) stop("k must be > 1 for this parameterization.", call. = FALSE)

  p <- .as_prob(p, log.p)
  if (!lower.tail) p <- 1 - p

  out <- rep.int(NA_real_, length(p))
  bad <- !is.finite(p) | p < 0 | p > 1
  out[bad] <- NaN
  out[p == 0] <- 0
  out[p == 1] <- Inf

  ok <- is.finite(p) & p > 0 & p < 1
  if (any(ok)) {
    out[ok] <- s * (((1 - p[ok])^(1 / (1 - k))) - 1)^(1 / c)
  }
  out
}

#' @rdname burr12
#' @export
rburr12 <- function(n, c, k, s) {
  .check_pos(c = c, s = s)
  if (any(k <= 1)) stop("k must be > 1 for this parameterization.", call. = FALSE)

  n <- as.integer(n)
  if (length(n) != 1 || n < 0) stop("n must be a nonnegative integer.", call. = FALSE)
  qburr12(stats::runif(n), c = c, k = k, s = s)
}

# -----------------------------
# Derived distribution (NEW) topic
# -----------------------------

#' Derived Burr-like distribution
#'
#' Density, distribution function, quantile function and random generation
#' for a derived Burr-like distribution with CDF
#' \deqn{F(x) = \left(1 + (T/x)^a \exp((T-x)/b)\right)^{-1},\quad x>0.}
#'
#' @param a,b,T Positive parameters.
#' @param tol Numerical tolerance passed to \code{\link[stats:uniroot]{uniroot}}.
#' @param maxiter Maximum iterations passed to \code{\link[stats:uniroot]{uniroot}}.
#' @inheritParams burr3
#' @name burrNEW
NULL

#' @rdname burrNEW
#' @export
dburrNEW <- function(x, a, b, T, log = FALSE) {
  .check_pos(a = a, b = b, T = T)
  x <- as.numeric(x)

  out <- rep.int(if (log) -Inf else 0, length(x))
  ok <- is.finite(x) & x > 0
  if (!any(ok)) return(out)

  xx <- x[ok]
  logg <- log(T) * a - a * log(xx) + (T - xx) / b
  g <- exp(logg)
  logdens <- log(xx / b + a) + logg - log(xx) - 2 * log1p(g)

  out[ok] <- if (log) logdens else exp(logdens)
  out
}

#' @rdname burrNEW
#' @export
pburrNEW <- function(q, a, b, T, lower.tail = TRUE, log.p = FALSE) {
  .check_pos(a = a, b = b, T = T)
  q <- as.numeric(q)

  p <- numeric(length(q))
  p[!is.finite(q)] <- NA_real_
  p[q <= 0] <- 0

  ok <- is.finite(q) & q > 0
  if (any(ok)) {
    xx <- q[ok]
    logg <- log(T) * a - a * log(xx) + (T - xx) / b
    g <- exp(logg)
    p[ok] <- 1 / (1 + g)
  }

  if (!lower.tail) p <- 1 - p
  if (log.p) p <- log(p)
  p
}

#' @rdname burrNEW
#' @export
qburrNEW <- function(p, a, b, T, lower.tail = TRUE, log.p = FALSE,
                     tol = 1e-10, maxiter = 200) {
  .check_pos(a = a, b = b, T = T)
  p <- .as_prob(p, log.p)
  if (!lower.tail) p <- 1 - p

  out <- rep.int(NA_real_, length(p))
  bad <- !is.finite(p) | p < 0 | p > 1
  out[bad] <- NaN
  out[p == 0] <- 0
  out[p == 1] <- Inf

  ok <- is.finite(p) & p > 0 & p < 1
  if (!any(ok)) return(out)

  for (i in which(ok)) {
    target <- p[i]
    f <- function(x) pburrNEW(x, a = a, b = b, T = T) - target

    lo <- .Machine$double.eps
    hi <- max(T, 1)

    fhi <- f(hi)
    it <- 0
    while (is.finite(fhi) && fhi < 0 && it < 60) {
      hi <- hi * 2
      fhi <- f(hi)
      it <- it + 1
    }

    if (!is.finite(fhi) || fhi < 0) {
      out[i] <- NA_real_
    } else {
      out[i] <- stats::uniroot(f, lower = lo, upper = hi,
                               tol = tol, maxiter = maxiter)$root
    }
  }
  out
}

#' @rdname burrNEW
#' @export
rburrNEW <- function(n, a, b, T) {
  .check_pos(a = a, b = b, T = T)
  n <- as.integer(n)
  if (length(n) != 1 || n < 0) stop("n must be a nonnegative integer.", call. = FALSE)
  qburrNEW(stats::runif(n), a = a, b = b, T = T)
}
