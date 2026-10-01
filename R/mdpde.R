
## -----------------------------------------------------------------------------
## (a) MDPDE for the linear regression model
##
##  Objective (up to an additive constant and the positive factor (2*pi)^(-a/2)):
##     f(b, s) = s^(-a) * [ (1+a)^(-1/2) - ((1+a)/a) * mean( exp(-a r_i^2/(2 s^2)) ) ]
##  with r_i = y_i - x_i' b.  At a = 0 it reduces to the Gaussian negative
##  log-likelihood, i.e. to least squares.
##
##  Two solvers are provided:
##   "irls"  (default) fixed-point / MM iteration, from the stationary equations
##              b   <- (X' W X)^{-1} X' W y,        W = diag(w_i),
##              s^2 <- mean(w r^2) / ( mean(w) - a (1+a)^(-3/2) ),
##              w_i  = exp(-a r_i^2 / (2 s^2)).
##           At a -> 0 the scale update returns mean(r^2), i.e. the MLE.
##   "optim" L-BFGS-B on (b, log s) with the analytic gradient supplied, box
##           constraints on log s.  Used as a fallback when IRLS stalls.
##
##  Initialisation: least squares for b, MAD of the LS residuals for s.  This is
##  a robust starting scale and is on the correct (unrescaled) footing.
##  Degeneracies: the scale is confined to [sigma_min, sigma_max]; the
##  denominator of the scale update is floored away from zero.  Both events are
##  counted and returned so that the driver can report how often they occur.
## -----------------------------------------------------------------------------
lmdpd <- function(y, X, alpha = 0.3,
                  method    = c("irls", "optim"),
                  tol       = 1e-8,
                  maxit     = 200L,
                  sigma_min = 1e-6,
                  sigma_max = 1e6,
                  init      = NULL) {

  method <- match.arg(method)
  X <- as.matrix(X); y <- as.numeric(y)
  n <- length(y); p <- ncol(X)

  ## --- initial values -------------------------------------------------------
  if (is.null(init)) {
    b0 <- tryCatch(qr.solve(X, y), error = function(e) rep(0, p))
    r0 <- y - X %*% b0
    s0 <- max(stats::mad(r0), 1e-3 * stats::sd(y), sigma_min)
  } else { b0 <- init$beta; s0 <- init$sigma }

  if (alpha == 0) {                       # least squares, closed form
    b <- tryCatch(qr.solve(X, y), error = function(e) rep(NA_real_, p))
    r <- y - X %*% b
    return(list(beta = as.numeric(b), sigma = sqrt(mean(r^2)),
                conv = 0L, iter = 0L, floored = FALSE))
  }

  cst <- alpha * (1 + alpha)^(-1.5)       # constant in the scale equation
  b <- b0; s <- s0; floored <- FALSE; it <- 0L; conv <- 1L

  if (method == "irls") {
    for (it in seq_len(maxit)) {
      r <- as.numeric(y - X %*% b)
      w <- exp(-alpha * r^2 / (2 * s^2))
      ## --- slope update: weighted least squares
      XW  <- X * w
      A   <- crossprod(XW, X)
      bn  <- tryCatch(solve(A, crossprod(XW, y)), error = function(e) NULL)
      if (is.null(bn)) { conv <- 52L; break }              # singular design
      ## --- scale update
      r   <- as.numeric(y - X %*% bn)
      w   <- exp(-alpha * r^2 / (2 * s^2))
      den <- mean(w) - cst
      if (den <= 1e-10) { den <- 1e-10; floored <- TRUE }  # implosion guard
      sn  <- sqrt(mean(w * r^2) / den)
      if (!is.finite(sn) || sn < sigma_min) { sn <- sigma_min; floored <- TRUE }
      if (sn > sigma_max) { sn <- sigma_max; floored <- TRUE }
      ## --- convergence on the relative change of (b, s)
      dd <- (sum(abs(bn - b)) + abs(sn - s)) / (sum(abs(b)) + abs(s) + 1e-10)
      b <- as.numeric(bn); s <- sn
      if (dd < tol) { conv <- 0L; break }
    }
    if (conv != 0L && conv != 52L) method <- "optim"       # fall through
  }

  if (method == "optim") {
    obj <- function(t) {
      bb <- t[seq_len(p)]; ss <- exp(t[p + 1])
      r  <- as.numeric(y - X %*% bb)
      ss^(-alpha) * ((1 + alpha)^(-0.5) -
                       ((1 + alpha) / alpha) * mean(exp(-alpha * r^2 / (2 * ss^2))))
    }
    gr <- function(t) {
      bb <- t[seq_len(p)]; ss <- exp(t[p + 1])
      r  <- as.numeric(y - X %*% bb)
      w  <- exp(-alpha * r^2 / (2 * ss^2))
      f  <- ss^(-alpha) * ((1 + alpha)^(-0.5) - ((1 + alpha) / alpha) * mean(w))
      gb <- -ss^(-alpha) * (1 + alpha) / ss^2 * colMeans(w * r * X)
      gu <- -alpha * f - ss^(-alpha) * (1 + alpha) / ss^2 * mean(w * r^2)
      c(gb, gu)
    }
    o <- try(stats::optim(c(b0, log(s0)), obj, gr, method = "L-BFGS-B",
                          lower = c(rep(-Inf, p), log(sigma_min)),
                          upper = c(rep( Inf, p), log(sigma_max)),
                          control = list(factr = 1e2, maxit = maxit)),
             silent = TRUE)
    if (inherits(o, "try-error"))
      return(list(beta = rep(NA_real_, p), sigma = NA_real_, conv = 99L,
                  iter = NA_integer_, floored = floored))
    b <- o$par[seq_len(p)]; s <- exp(o$par[p + 1]); conv <- o$convergence
    it <- o$counts[1]
  }

  list(beta = as.numeric(b), sigma = s, conv = as.integer(conv),
       iter = as.integer(it), floored = floored)
}



## -----------------------------------------------------------------------------
## (b) MDPDE for the intercept-only LMM   y_i = mu 1 + Z_i b_i + e_i
##
##  Cluster-wise density power divergence for independent, non-homogeneous
##  multivariate normals (Saraceno et al., 2024, in spirit):
##     V_a(y_i) = (2 pi)^(-a n_i/2) |V_i|^(-a/2) *
##                [ (1+a)^(-n_i/2) - (1 + 1/a) exp(-a/2 * r_i' V_i^{-1} r_i) ]
##  with V_i = sigma^2 (I + Z_i P Z_i'), P = Psi / sigma^2.
##  Parameterisation: mu, log sigma^2 and the log-Cholesky factor of Psi, so the
##  optimisation is unconstrained and Psi stays positive definite.
##  Returns P = Psi / sigma^2, which is exactly the proxy matrix required.
##  This is a single low-dimensional fit, so its cost does not scale with p.
## -----------------------------------------------------------------------------
mdpde_lmm_icept <- function(y, Z, ID, alpha = 0.3, maxit = 500L) {
  Z <- as.matrix(Z); q <- ncol(Z)
  idx <- split(seq_along(y), ID)
  nl  <- q * (q + 1) / 2

  unpack <- function(t) {
    mu <- t[1]; s2 <- exp(t[2])
    L  <- matrix(0, q, q); L[lower.tri(L, diag = TRUE)] <- t[3:(2 + nl)]
    diag(L) <- exp(diag(L))
    list(mu = mu, s2 = s2, Psi = tcrossprod(L))
  }
  obj <- function(t) {
    par <- unpack(t); s2 <- par$s2; Psi <- par$Psi
    tot <- 0
    for (g in idx) {
      Zi <- Z[g, , drop = FALSE]; ni <- length(g)
      Vi <- s2 * diag(ni) + Zi %*% Psi %*% t(Zi)
      ch <- tryCatch(chol(Vi), error = function(e) NULL)
      if (is.null(ch)) return(1e10)
      ldet <- 2 * sum(log(diag(ch)))
      r    <- y[g] - par$mu
      quad <- sum(backsolve(ch, r, transpose = TRUE)^2)
      lk   <- -0.5 * alpha * (ni * log(2 * pi) + ldet)      # log |V|^(-a/2) etc.
      tot  <- tot + exp(lk) * ((1 + alpha)^(-ni / 2) -
                                 (1 + 1 / alpha) * exp(-0.5 * alpha * quad))
    }
    tot / length(idx)
  }
  t0 <- c(mean(y), log(max(stats::var(y) / 2, 1e-3)),
          c(log(rep(sqrt(max(stats::var(y) / 4, 1e-3)), q)), rep(0, nl - q))[
            order(c(seq_len(q), rep(q + 1, nl - q)))][seq_len(nl)])
  ## simpler, safer starting value for the log-Cholesky block:
  L0 <- diag(q) * sqrt(max(stats::var(y) / 4, 1e-3))
  diag(L0) <- log(diag(L0))
  t0 <- c(mean(y), log(max(stats::var(y) / 2, 1e-3)),
          L0[lower.tri(L0, diag = TRUE)])
  o <- try(stats::optim(t0, obj, method = "Nelder-Mead",
                        control = list(maxit = maxit, reltol = 1e-9)),
           silent = TRUE)
  if (inherits(o, "try-error"))
    return(list(P = NULL, sigma2 = NA_real_, Psi = NULL, conv = 99L))
  par <- unpack(o$par)
  list(P = par$Psi / par$s2, sigma2 = par$s2, Psi = par$Psi,
       conv = as.integer(o$convergence))
}
