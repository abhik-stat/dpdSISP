## =============================================================================
##  Compute proxy matrices, whitening, and whitening-error diagnostics
## =============================================================================


## -----------------------------------------------------------------------------
## Inverse symmetric square root of a p.d. matrix, with a floor on the
## eigenvalues (this is the constant C_m^tilde of assumption (P0)).
## -----------------------------------------------------------------------------
inv_sqrt_sym <- function(V, eps = 1e-10) {
  e <- eigen((V + t(V)) / 2, symmetric = TRUE)
  d <- pmax(e$values, eps)
  e$vectors %*% (t(e$vectors) / sqrt(d))
}


## -----------------------------------------------------------------------------
## Block-diagonal whitening.  Returns the transformed response/design and,
## optionally, the realised whitening error against a known truth.
##
##  P     q x q proxy matrix (the estimate of Psi / sigma^2)
##  Ptrue q x q true value, if available (oracle diagnostics only)
## -----------------------------------------------------------------------------
whiten_by_cluster <- function(y, X, Z, ID, P, Ptrue = NULL) {
  X <- as.matrix(X); Z <- as.matrix(Z)
  ya <- numeric(length(y)); Xa <- matrix(0, nrow(X), ncol(X))
  err_S <- 0; err_V <- 0
  for (g in split(seq_along(y), ID)) {
    Zi  <- Z[g, , drop = FALSE]; ni <- length(g)
    Si  <- Zi %*% P %*% t(Zi) + diag(ni)
    iS  <- inv_sqrt_sym(Si)
    ya[g]  <- iS %*% y[g]
    Xa[g, ] <- iS %*% X[g, , drop = FALSE]
    if (!is.null(Ptrue)) {
      St <- Zi %*% Ptrue %*% t(Zi) + diag(ni)
      err_S <- max(err_S, norm(Si - St, "2"))              # ||Sigma_hat - Sigma||
      Vi <- iS %*% St %*% iS
      err_V <- max(err_V, norm(Vi - diag(ni), "2"))        # ||V_i - I||
    }
  }
  list(y = ya, X = Xa, err_Sigma = err_S, err_V = err_V)
}


## -----------------------------------------------------------------------------
## D_Z = diag(root-mean-square of each column of Z).
## Used by: (i) standardize_Z(), just below, and (ii) the "logn" benchmark
## proxy further down (Fan & Li's 2012 fixed choice).
##
## -----------------------------------------------------------------------------
DZ_scale <- function(Z) {
  s <- sqrt(colMeans(as.matrix(Z)^2))
  s[s < 1e-8] <- 1
  s
}


## -----------------------------------------------------------------------------
## Z-standardisation: Z_i <- Z_i D_Z^{-1} (paper Section 2), applied to
## whatever raw Z a pipeline has just obtained
##
## Idempotent: applying this to an already-standardised Z (root-mean-square
## already 1 on every non-intercept column) leaves it unchanged, so calling
## it more than once on the same Z along a code path is harmless if it ever
## happens, just redundant.
## -----------------------------------------------------------------------------
standardize_Z <- function(Z) {
  Z <- as.matrix(Z)
  n <- nrow(Z); q <- ncol(Z)
  has_icept <- isTRUE(all.equal(Z[, 1], rep(1, n)))
  z_cols <- if (has_icept) { if (q > 1L) 2:q else integer(0) } else seq_len(q)
  if (length(z_cols)) {
    DZ <- DZ_scale(Z[, z_cols, drop = FALSE])
    Z[, z_cols] <- sweep(Z[, z_cols, drop = FALSE], 2, DZ, "/")
  }
  Z
}


## -----------------------------------------------------------------------------
## make_proxy(): returns a q x q matrix P.
##
##  type = "oracle"   P = Psi_true / sigma2_true        (delta_n = 0)
##         "cvP"      P = a_opt * I_q, a_opt by K-fold CV (paper's cv-P)
##         "I0"       P = Psi_hat / sigma2_hat, intercept-only REML
##         "robP"     P from the intercept-only MDPDE (robust version of I0)
##         "none"     P = 0  (no whitening; the P = O_q benchmark)
##         "scaled"   P = ctrl$fac * oracle   (deliberate misspecification)
##         "logn"     P = log(n) * D_Z^{-2}   (the fixed choice of Fan & Li 2012)
##
##  idx  : optional subset of rows on which the proxy is to be estimated.  Pass
##         the proxy half of a cluster split here to make Sigma_hat independent
##         of the screening statistics.
## -----------------------------------------------------------------------------
make_proxy <- function(type, y, X, Z, ID, idx = NULL,
                       truth = NULL, alpha = 0.3, ctrl = list()) {
  Z <- as.matrix(Z); q <- ncol(Z)
  if (is.null(idx)) idx <- seq_along(y)
  yy <- y[idx]; ZZ <- Z[idx, , drop = FALSE]; II <- ID[idx]
  # Dz <- DZ_scale(Z)

  switch(type,
    "oracle" = {
      if (is.null(truth)) stop("oracle proxy needs truth = list(Psi=, sigma2=)")
      truth$Psi / truth$sigma2
    },
    "none"   = matrix(0, q, q),
    "logn"   = log(length(y)) * diag(q),
    "scaled" = {
      if (is.null(truth)) stop("scaled proxy needs truth")
      (ctrl$fac %||% 4) * truth$Psi / truth$sigma2
    },
    "cvP"    = cvP_proxy(yy, X[idx, , drop = FALSE], ZZ, II, q,
                         a_grid = ctrl$a_grid %||% seq(0, 20, by = 1),
                         nfold  = ctrl$nfold  %||% 5),
    "I0"     = I0_proxy(yy, ZZ, II, q),
    "robP"   = {
      f <- mdpde_lmm_icept(yy, ZZ, II, alpha = alpha)
      if (is.null(f$P)) {
        I0_proxy(yy, ZZ, II, q)
      } else {
        ev <- eigen(f$P, symmetric = TRUE, only.values = TRUE)$values
        cond_max <- ctrl$robP_cond_max %||% 200
        if (min(ev) <= 1e-6 || max(ev) / max(min(ev), 1e-10) > cond_max)
          I0_proxy(yy, ZZ, II, q)
        else f$P
      }
    },
    stop("unknown proxy type: ", type)
  )
}


## --- intercept-only REML proxy (I0-P) ----------------------------------------
I0_proxy <- function(y, Z, ID, q) {
  Z <- as.matrix(Z)
  if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(q) - 1L)
  dat <- data.frame(as.data.frame(Z), y = as.numeric(y), ID = as.factor(ID))
  fml <- if (q > 1)
    stats::as.formula(paste0("y ~ (1 + ",
                             paste0(colnames(Z)[-1], collapse = "+"), " | ID)"))
  else stats::as.formula("y ~ (1 | ID)")
  fit <- try(lme4::lmer(fml, data = dat, REML = TRUE), silent = TRUE)
  if (inherits(fit, "try-error")) return(diag(0, q))
  Psi <- as.matrix(lme4::VarCorr(fit)$ID)[seq_len(q), seq_len(q), drop = FALSE]
  Psi / (lme4::getME(fit, "sigma")^2)
}


## --- cross-validated scalar proxy (cv-P) --------------------------------------
##  P_cv = a_opt * I_q (paper's Eq. for cv-P), a_opt minimising the K-fold
##  prediction error of a lasso fit on the whitened data. Folds are taken over
##  CLUSTERS, not records, so that the whitening of a test cluster never uses
##  training information from the same cluster; and the whitening is block
##  diagonal. P_cv = a_opt * I_q is used directly, with no further rescaling here.
##  
cvP_proxy <- function(y, X, Z, ID, q, a_grid, nfold) {
  Z <- as.matrix(Z); X <- as.matrix(X)
  clus  <- unique(ID)
  folds <- split(clus, sample(rep(seq_len(nfold), length.out = length(clus))))
  err   <- numeric(length(a_grid))
  for (i in seq_along(a_grid)) {
    P <- a_grid[i] * diag(q)
    e <- numeric(length(folds))
    for (k in seq_along(folds)) {
      te <- which(ID %in% folds[[k]]); tr <- setdiff(seq_along(y), te)
      wtr <- whiten_by_cluster(y[tr], X[tr, , drop = FALSE],
                               Z[tr, , drop = FALSE], ID[tr], P)
      wte <- whiten_by_cluster(y[te], X[te, , drop = FALSE],
                               Z[te, , drop = FALSE], ID[te], P)
      fit <- try(glmnet::cv.glmnet(wtr$X, wtr$y, alpha = 1, standardize = FALSE),
                 silent = TRUE)
      if (inherits(fit, "try-error")) { e[k] <- NA; next }
      pr  <- as.numeric(predict(fit, newx = wte$X, s = "lambda.min"))
      e[k] <- mean((wte$y - pr)^2)          # error on the whitened scale
    }
    err[i] <- mean(e, na.rm = TRUE)
  }
  a_opt <- a_grid[which.min(err)]
  attr_P <- a_opt * diag(q)
  attr(attr_P, "a_opt") <- a_opt
  attr(attr_P, "cv_err") <- err
  attr_P
}


## -----------------------------------------------------------------------------
## Diagnostic reported for every configuration: the realised
## spectral accuracy of the proxy, which is the quantity assumption (P0_delta)
## constrains.  Requires the truth, so it is computed in simulations only.
## -----------------------------------------------------------------------------
proxy_diagnostics <- function(Z, ID, P, truth) {
  Ptrue <- truth$Psi / truth$sigma2
  eS <- 0; eV <- 0
  for (g in split(seq_along(ID), ID)) {
    Zi <- as.matrix(Z)[g, , drop = FALSE]; ni <- length(g)
    Sh <- Zi %*% P     %*% t(Zi) + diag(ni)
    St <- Zi %*% Ptrue %*% t(Zi) + diag(ni)
    iS <- inv_sqrt_sym(Sh)
    eS <- max(eS, norm(Sh - St, "2"))
    eV <- max(eV, norm(iS %*% St %*% iS - diag(ni), "2"))
  }
  c(err_Sigma = eS, err_V = eV)
}
