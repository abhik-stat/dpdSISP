## =============================================================================
##  0screening procedures
##
##   screen_dpd()    DPD-SISP / DPD-CSISP with any proxy, scale-invariant rank
##   screen_bench()  benchmark marginal LMM screening (MLE, REML, HGD(gamma))
##   screen_isis()   DPD-ISISP, iterative (bug-fixed, see notes)
##   screen_tpcc()   thresholded partial correlation, conditional (TPCc)
##
##  CONVENTION.  X has the intercept in column 1; screening is over columns
##  2..p.  All functions return indices IN X-COLUMN SPACE (i.e. values in
##  2..p), ranked by decreasing importance.  This differs from the previous
##  code, which returned 1..(p-1) and required a "+1" at the call site; doing it
##  once here removes an easy off-by-one.
##
## -----------------------------------------------------------------------------
screen_dpd <- function(y, X, Z, ID, alpha = 0.3, proxy = "cvP",
                       Cond = NULL, truth = NULL, weighted = FALSE,
                       split = FALSE, ctrl = list(), diagnostics = TRUE) {

  X <- as.matrix(X); Z <- as.matrix(Z)
  n <- length(y); p <- ncol(X)
  scr_cols  <- setdiff(2:p, Cond)          # columns to be screened
  cond_cols <- c(1L, Cond)                 # intercept + conditioning set

  ## --- proxy, optionally on an independent half of the clusters -------------
  if (split && proxy %in% c("cvP", "I0", "robP")) {
    sp   <- split_clusters(ID)
    Pmat <- make_proxy(proxy, y, X, Z, ID, idx = sp$prox,
                       truth = truth, alpha = alpha, ctrl = ctrl)
    keep <- sp$screen                      # screen on the other half
  } else {
    Pmat <- make_proxy(proxy, y, X, Z, ID, truth = truth,
                       alpha = alpha, ctrl = ctrl)
    keep <- seq_len(n)
  }

  ## --- whitening (block diagonal) -------------------------------------------
  wh <- whiten_by_cluster(y[keep], X[keep, , drop = FALSE],
                          Z[keep, , drop = FALSE], ID[keep], Pmat,
                          Ptrue = if (diagnostics && !is.null(truth))
                            truth$Psi / truth$sigma2 else NULL)
  ya <- wh$y; Xa <- wh$X

  ## --- centre and scale the whitened covariates ------------------------------
  ctr <- colMeans(Xa[, scr_cols, drop = FALSE])
  sdv <- apply(Xa[, scr_cols, drop = FALSE], 2, stats::sd)
  sdv[sdv < 1e-10] <- 1

  bet <- rep(NA_real_, length(scr_cols))
  nconv <- 0L; nfloor <- 0L
  fitter <- if (weighted) lmdpd_weighted else lmdpd
  for (k in seq_along(scr_cols)) {
    xj  <- (Xa[, scr_cols[k]] - ctr[k]) / sdv[k]        # x_tilde_j_std, Eq. (7)
    Xk  <- cbind(Xa[, cond_cols, drop = FALSE], xj)
    fit <- try(fitter(ya, Xk, alpha), silent = TRUE)    # Eq. (8)/(9): MDPDE on
                                                          # the standardised fit
    if (inherits(fit, "try-error") || fit$conv != 0L) nconv <- nconv + 1L
    if (!inherits(fit, "try-error")) {
      bet[k] <- fit$beta[length(fit$beta)]      # beta_hat_j1_alpha, Eq. (8)-(9)
      if (isTRUE(fit$floored)) nfloor <- nfloor + 1L
    }
  }

  ## beta_hat is the paper's own |beta_hat_j1_alpha| -- the screening statistic
  ## itself (Eq. (11)-(12)), not a separately-defined quantity. raw is kept
  ## only as a diagnostic: the slope a fit on the Unstandardised covariate
  ## would have given (i.e. beta_hat / sd), which is scale-DEPENDENT and is
  ## not what the paper's procedure ranks by; it is reported so a check can
  ## confirm the two rankings coincide when covariates already share a common
  ## scale, and diverge otherwise.
  beta_hat <- abs(bet)
  raw      <- abs(bet) / sdv
  ord      <- order(beta_hat, decreasing = TRUE, na.last = NA)

  list(order   = scr_cols[ord],
       beta_hat = beta_hat[ord],
       U       = beta_hat[ord],            # alias, kept for backward compatibility
       order_raw = scr_cols[order(raw, decreasing = TRUE, na.last = NA)],
       stat    = setNames(beta_hat, scr_cols),
       proxy   = Pmat,
       a_opt   = attr(Pmat, "a_opt") %||% NA_real_,
       err_Sigma = wh$err_Sigma, err_V = wh$err_V,
       n_nonconv = nconv, n_floored = nfloor)
}


screen_dpd_lm <- function(y, X, alpha = 0.3, 
                       Cond = NULL, truth = NULL, weighted = FALSE,
                       split = TRUE, ctrl = list(), diagnostics = TRUE) {
  
  X <- as.matrix(X)
  n <- length(y); p <- ncol(X)
  scr_cols  <- setdiff(2:p, Cond)          # columns to be screened
  cond_cols <- c(1L, Cond)                 # intercept + conditioning set
  

  ## --- centre and scale the whitened covariates ------------------------------
  ctr <- colMeans(X[, scr_cols, drop = FALSE])
  sdv <- apply(X[, scr_cols, drop = FALSE], 2, stats::sd)
  sdv[sdv < 1e-10] <- 1
  
  bet <- rep(NA_real_, length(scr_cols))
  nconv <- 0L; nfloor <- 0L
  fitter <- if (weighted) lmdpd_weighted else lmdpd
  for (k in seq_along(scr_cols)) {
    xj  <- (X[, scr_cols[k]] - ctr[k]) / sdv[k]        # x_tilde_j_std, Eq. (7)
    Xk  <- cbind(X[, cond_cols, drop = FALSE], xj)
    fit <- try(fitter(y, Xk, alpha), silent = TRUE)    # Eq. (8)/(9): MDPDE on
    # the standardised fit
    if (inherits(fit, "try-error") || fit$conv != 0L) nconv <- nconv + 1L
    if (!inherits(fit, "try-error")) {
      bet[k] <- fit$beta[length(fit$beta)]      # beta_hat_j1_alpha, Eq. (8)-(9)
      if (isTRUE(fit$floored)) nfloor <- nfloor + 1L
    }
  }
  
  ## beta_hat is the paper's own |beta_hat_j1_alpha| -- the screening statistic
  ## itself (Eq. (11)-(12)), not a separately-defined quantity. raw is kept
  ## only as a diagnostic: the slope a fit on the UNstandardised covariate
  ## would have given (i.e. beta_hat / sd), which is scale-DEPENDENT and is
  ## not what the paper's procedure ranks by; it is reported so a check can
  ## confirm the two rankings coincide when covariates already share a common
  ## scale, and diverge otherwise.
  beta_hat <- abs(bet)
  raw      <- abs(bet) / sdv
  ord      <- order(beta_hat, decreasing = TRUE, na.last = NA)
  
  list(order   = scr_cols[ord],
       beta_hat = beta_hat[ord],
       U       = beta_hat[ord],            # alias, kept for backward compatibility
       order_raw = scr_cols[order(raw, decreasing = TRUE, na.last = NA)],
       stat    = setNames(beta_hat, scr_cols),
       proxy   = NA_real_,
       a_opt   = NA_real_,
       err_Sigma = NA_real_, err_V = NA_real_,
       n_nonconv = nconv, n_floored = nfloor)
}


## -----------------------------------------------------------------------------
## Benchmark screening: marginal mixed models fitted by ML, REML or HGD(gamma).
## Uses LMM_MLE / LMM_REML / LMM_HGD.
## -----------------------------------------------------------------------------
screen_bench <- function(y, X, Z, ID, method = c("MLE", "REML", "HGD"),
                         gamma = 0.3, Cond = NULL, maxItr = 100L) {
  method <- match.arg(method)
  X <- as.matrix(X); p <- ncol(X)
  scr_cols  <- setdiff(2:p, Cond)
  cond_cols <- c(1L, Cond)
  f <- switch(method, MLE = LMM_MLE, REML = LMM_REML, HGD = LMM_HGD)
  est <- rep(NA_real_, length(scr_cols))
  for (k in seq_along(scr_cols)) {
    Xk  <- X[, c(cond_cols, scr_cols[k]), drop = FALSE]
    fit <- try(f(as.numeric(y), Xk, Z, ID, gam = gamma, maxitr = maxItr),
               silent = TRUE)
    if (!inherits(fit, "try-error")) est[k] <- fit$Beta[length(cond_cols) + 1L]
  }
  ## same standardization as above, on the non-transformed covariates
  sdv <- apply(X[, scr_cols, drop = FALSE], 2, stats::sd); sdv[sdv < 1e-10] <- 1
  U   <- abs(est) * sdv
  ord <- order(U, decreasing = TRUE, na.last = NA)
  list(order = scr_cols[ord], U = U[ord], stat = setNames(U, scr_cols))
}


## -----------------------------------------------------------------------------
## Iterative DPD-SISP (DPD-ISISP).
##
## Fixes relative to DPD_ISISP_LMM.r:
##   (i)   batch[itr] was used before itr existed;
##   (ii)  `break` was called outside any loop (a hard error in R);
##   (iii) the residual step used the full parameter vector, whose last element
##         is the scale, in  ya - Xa[, IR2] %*% est10  -- a dimension mismatch;
##   (iv)  whitening is now block diagonal;
##   (v)   the marginal statistic is |beta_hat_j1_alpha|, fit directly on the
##         standardised residualised covariate (paper Eq. (7)-(8), applied at
##         each iteration to the current residual in place of Y).
## -----------------------------------------------------------------------------
screen_isis <- function(y, X, Z, ID, d, alpha = 0.3, proxy = "cvP",
                        Cond = NULL, truth = NULL, maxItr = 5L,
                        batch = NULL, ctrl = list()) {

  X <- as.matrix(X); p <- ncol(X); n <- length(y)
  if (is.null(batch)) batch <- max(1L, ceiling(d / maxItr))
  Pmat <- make_proxy(proxy, y, X, Z, ID, truth = truth, alpha = alpha, ctrl = ctrl)
  wh <- whiten_by_cluster(y, X, Z, ID, Pmat)
  ya <- wh$y; Xa <- wh$X

  std <- function(cols) {
    M <- Xa[, cols, drop = FALSE]
    s <- apply(M, 2, stats::sd); s[s < 1e-10] <- 1
    sweep(sweep(M, 2, colMeans(M)), 2, s, "/")
  }
  ## fits the marginal MDPDE directly on the standardised covariate (Eq. (7)),
  ## so f$beta[2] is already beta_hat_j1_alpha -- no rescaling step needed.
  marg <- function(resp, cols) {
    S <- std(cols); out <- rep(NA_real_, length(cols))
    for (k in seq_along(cols)) {
      f <- try(lmdpd(resp, cbind(1, S[, k]), alpha), silent = TRUE)
      if (!inherits(f, "try-error")) out[k] <- abs(f$beta[2])
    }
    out
  }

  active <- c(1L, Cond)                         # intercept + conditioning set
  cand   <- setdiff(2:p, Cond)
  sel    <- integer(0)

  u  <- marg(ya, cand)
  ad <- cand[order(u, decreasing = TRUE, na.last = NA)][seq_len(min(batch, sum(!is.na(u))))]
  sel <- ad; itr <- 1L

  while (length(sel) < d && itr < maxItr) {
    itr <- itr + 1L
    cols <- c(active, sel)
    fit  <- try(lmdpd(ya, Xa[, cols, drop = FALSE], alpha), silent = TRUE)
    if (inherits(fit, "try-error") || anyNA(fit$beta)) break
    res  <- ya - as.numeric(Xa[, cols, drop = FALSE] %*% fit$beta)   # slopes only
    cand <- setdiff(cand, sel)
    if (!length(cand)) break
    u  <- marg(res, cand)
    if (all(is.na(u))) break
    k  <- min(batch, sum(!is.na(u)), d - length(sel))
    ad <- cand[order(u, decreasing = TRUE, na.last = NA)][seq_len(k)]
    if (!length(ad)) break
    sel <- c(sel, ad)
  }
  list(order = sel[seq_len(min(length(sel), d))], iter = itr)
}


## =============================================================================
##  TPCc -- conditional thresholded partial correlation
##  Alabiso & Shang (2023), Comm. Statist. Theory Methods 52, 6355-6380,
##  which applies the TPC algorithm of Li, Liu & Lou (2017), Statist. Sinica
##  27(3), 983-996, to the partial correlation between covariates and response
##  CONDITIONAL ON THE RANDOM EFFECTS.
##
##  The TPC part below follows Algorithm 1 and equations (2.6)-(2.9) of
##  Li, Liu & Lou (2017) exactly:
##
##    kappa_hat = (1/p) sum_j [ m4_j / (3 m2_j^2) - 1 ]                    (2.8)
##    T(a,n,k,|S|) = (e^{2w} - 1)/(e^{2w} + 1),
##        w = sqrt(1+k) * Phi^{-1}(1 - a/2) / sqrt(n - 1 - |S|)            (2.7)
##    A[1]   = { j : |rho_hat(y, x_j)| > T(a,n,k,0) }
##    A[m]   = { j in A[m-1] : |rho_hat(y, x_j | x_S)| > T(a,n,k,|S|)
##                             for all S subset A[m-1]\{j}, |S| = m-1 }
##    stop at m_reach = min{ m : |A[m]| <= m }
##
##  For kappa = 0 this is the PC-simple algorithm of Buhlmann, Kalisch and
##  Maathuis (2010), as it should be.
##
##  The CONDITIONING STEP is our reading of "partial correlation conditional 
##  on the random effects":
##     condition = "blup"   subtract the BLUPs of an intercept-only LMM from the
##                          response, y_adj = y - Z_i bhat_i, then run TPC on
##                          (y_adj, X).  This is the literal reading and the
##                          default.
##     condition = "whiten" run TPC on the proxy-whitened data, which is the
##                          construction used elsewhere in this paper.
##     condition = "none"   plain TPC, ignoring the random effects (the TPC that
##                          Alabiso & Shang report TPCc to beat).
##
##  TWO PRACTICAL CAPS, both reported in the return value:
##   max_set   the first-step set is truncated to its strongest `max_set`
##             members.  Step m of the exact algorithm examines
##             choose(|A[m-1]|-1, m-1) subsets per variable, which is infeasible
##             when the first step keeps hundreds of covariates -- as it does at
##             p = 20000 under contamination.  Li et al. never face this because
##             their |A[1]| is small.
##   max_order the maximum conditioning order m.
##  Both caps only ever make TPCc keep MORE variables, never fewer, so they
##  cannot manufacture an advantage for our method.
##
##  RANKING.  TPC returns a set, not an ordering, whereas TPR / MinMS need a
##  ranking.  We rank by the step at which a covariate was eliminated (survivors
##  first), breaking ties by |marginal correlation|.  This is stated in the
##  paper alongside the comparison.
## =============================================================================

## excess-kurtosis estimate, eq. (2.8) of Li, Liu & Lou (2017) 
tpc_kurtosis <- function(X) {
  k <- apply(X, 2, function(x) {
    x <- x - mean(x); m2 <- mean(x^2); m4 <- mean(x^4)
    if (m2 <= 0) NA_real_ else m4 / (3 * m2^2) - 1
  })
  max(mean(k, na.rm = TRUE), -0.99)      # (D3) requires kappa > -1
}

## rejection threshold, eq. (3.7) of Alabiso & Shang (2023) 
tpc_threshold <- function(alpha, n, kappa, s) {
  df <- max(n - 1 - s, 1)
  w  <- sqrt(1 + kappa) * stats::qnorm(1 - alpha / 2) / sqrt(df)
  (exp(2 * w) - 1) / (exp(2 * w) + 1)
}

## partial correlation via residualization, Definition 3.1 
## XS may be NULL/zero-column (unconditional correlation).
pcor_given <- function(y, Xc, XS) {
  if (is.null(XS) || !ncol(XS)) {
    ry <- y - mean(y); rX <- sweep(Xc, 2, colMeans(Xc))
  } else {
    q  <- qr(cbind(1, XS))
    ry <- qr.resid(q, y); rX <- qr.resid(q, Xc)
  }
  sy <- sqrt(sum(ry^2)); sx <- sqrt(colSums(rX^2))
  ok <- sx > 1e-10 & sy > 1e-10
  out <- rep(0, ncol(Xc))
  out[ok] <- as.numeric(crossprod(rX[, ok, drop = FALSE], ry)) / (sx[ok] * sy)
  out
}

## -----------------------------------------------------------------------------
## Core engine: Algorithm 1, Steps 1-3, on a given response `resp` against
## candidate columns `Xcand`, with `Xforce` (possibly NULL) always included
## in every conditioning set tested. Used for BOTH Step 0a (on the raw
## response) and the final Steps 1-3 (on y_b) -- same algorithm, two calls.
## `max_set`/`max_order` are our own tractability caps (see header of the
## original screen_tpcc block above in the file for the full rationale);
## they can only enlarge the returned active set, never shrink it.
## -----------------------------------------------------------------------------
run_tpc_steps <- function(resp, Xcand, Xforce, alpha_level, max_set, max_order) {
  resp  <- as.numeric(scale(resp))
  Xcand <- scale(Xcand); Xcand[!is.finite(Xcand)] <- 0
  if (!is.null(Xforce) && ncol(Xforce)) {
    Xforce <- scale(Xforce); Xforce[!is.finite(Xforce)] <- 0
  } else Xforce <- NULL
  n <- length(resp)
  n_force <- if (is.null(Xforce)) 0L else ncol(Xforce)
  kappa <- tpc_kurtosis(Xcand)
  
  ## Step 1 (conditioned on Xforce alone, if supplied)
  rho1 <- abs(pcor_given(resp, Xcand, Xforce))
  T0   <- tpc_threshold(alpha_level, n, kappa, n_force)
  A    <- which(rho1 > T0)
  truncated <- FALSE
  if (length(A) > max_set) {
    A <- A[order(rho1[A], decreasing = TRUE)][seq_len(max_set)]
    truncated <- TRUE
  }
  drop_step <- rep(NA_integer_, ncol(Xcand))
  drop_step[setdiff(seq_len(ncol(Xcand)), A)] <- 1L
  
  ## Steps 2, 3, ... (conditioned on Xforce plus a growing subset of size m-1)
  m <- 1L
  while (length(A) > m && m < max_order) {
    m <- m + 1L
    Tm <- tpc_threshold(alpha_level, n, kappa, n_force + m - 1L)
    alive <- A
    for (S in utils::combn(A, m - 1L, simplify = FALSE)) {
      test <- setdiff(alive, S)
      if (!length(test)) next
      Sfull <- cbind(Xcand[, S, drop = FALSE], Xforce)
      r <- abs(pcor_given(resp, Xcand[, test, drop = FALSE], Sfull))
      kill <- test[r <= Tm]
      if (length(kill)) {
        drop_step[kill] <- m
        alive <- setdiff(alive, kill)
      }
      if (length(alive) <= m - 1L) break
    }
    if (identical(sort(alive), sort(A))) { A <- alive; break }   # m_reach
    A <- alive
  }
  
  list(active = A, drop_step = drop_step, rho1 = rho1, kappa = kappa,
       threshold = T0, m_reach = m, truncated = truncated)
}


## -----------------------------------------------------------------------------
## screen_tpcc(): Alabiso & Shang (2023)'s TPCc, the two-stage version.
##
## `condition = "blup"` is now the paper's actual mechanism (Step 0a/0b/0c
## below, matching Section 3.1 exactly); `condition = "whiten"` and
## `"none"` remain our own alternatives, NOT described in the paper
## -----------------------------------------------------------------------------
screen_tpcc <- function(y, X, Z, ID, d = NULL, Cond = NULL,
                        condition = c("blup", "whiten", "none"),
                        alpha_level = 0.05, proxy = "I0", truth = NULL,
                        max_set = 60L, max_order = 3L) {
  
  condition <- match.arg(condition)
  X <- as.matrix(X); Z <- as.matrix(Z); p <- ncol(X)
  cand   <- setdiff(2:p, Cond)
  Xforce <- if (is.null(Cond)) NULL else X[, Cond, drop = FALSE]
  init_active <- NULL
  
  if (condition == "blup") {
    ## ---- Step 0a: initial active set via (unconditional) TPC on the RAW
    ## response, full candidate set (paper's Sec. 3.1: "we use the TPC
    ## algorithm on the unconditioned response to estimate an initial
    ## active set").
    step0 <- run_tpc_steps(y, X[, cand, drop = FALSE], Xforce,
                           alpha_level, max_set, max_order)
    init_active <- cand[step0$active]
    
    ## ---- Step 0b: REML fit of y ~ Cond + X[,init_active] + (Z|ID) -- WITH
    ## the initial active set (and Cond) as fixed effects, matching the
    ## paper's eq. (3.1) exactly. This is the step that was wrong before:
    ## the previous version fit an intercept-only random-effects model here.
    ## ---- Step 0c: BLUP from that model; y_b = y - Z*bhat (random part only).
    fe_cols <- c(Cond, init_active)
    q <- ncol(Z)
    if (is.null(colnames(Z))) colnames(Z) <- paste0("Z", seq_len(q) - 1L)
    
    if (!length(fe_cols)) {
      ## degenerate case: nothing survived Step 0a and no Cond supplied --
      ## fall back to an intercept-only random-effects model, since there is
      ## no fixed-effect information to include (this is the ONLY situation
      ## in which the previous behavior was actually correct, by necessity).
      dat <- data.frame(as.data.frame(Z), y = as.numeric(y), ID = as.factor(ID))
      fml <- if (q > 1)
        stats::as.formula(paste0("y ~ (1 + ",
                                 paste0(colnames(Z)[-1], collapse = "+"), " | ID)"))
      else stats::as.formula("y ~ (1 | ID)")
    } else {
      Xfe <- as.data.frame(X[, fe_cols, drop = FALSE])
      names(Xfe) <- paste0("fe", seq_along(fe_cols))
      dat <- data.frame(as.data.frame(Z), Xfe, y = as.numeric(y), ID = as.factor(ID))
      fe_rhs <- paste(names(Xfe), collapse = " + ")
      re_rhs <- if (q > 1)
        paste0("(1 + ", paste0(colnames(Z)[-1], collapse = "+"), " | ID)")
      else "(1 | ID)"
      fml <- stats::as.formula(paste("y ~", fe_rhs, "+", re_rhs))
    }
    
    fit <- try(lme4::lmer(fml, data = dat, REML = TRUE), silent = TRUE)
    if (inherits(fit, "try-error")) {
      yy <- y   # degrade to unconditioned rather than fail outright
    } else {
      B   <- as.matrix(lme4::ranef(fit)$ID)          # m x q BLUPs
      idx <- match(ID, sort(unique(ID)))
      yy  <- y - rowSums(Z * B[idx, , drop = FALSE])  # random part ONLY
    }
    
  } else if (condition == "whiten") {
    ## NOT the paper's mechanism -- our own alternative, kept for comparison.
    P  <- make_proxy(proxy, y, X, Z, ID, truth = truth)
    wh <- whiten_by_cluster(y, X, Z, ID, P)
    yy <- wh$y
    X  <- wh$X                                        # whitening transforms X too
    Xforce <- if (is.null(Cond)) NULL else X[, Cond, drop = FALSE]
    
  } else { yy <- y }                                  # condition == "none"
  
  ## ---- Steps 1-3 proper: TPCc on yb (or plain TPC on y if condition="none")
  final <- run_tpc_steps(yy, X[, cand, drop = FALSE], Xforce,
                         alpha_level, max_set, max_order)
  
  key <- ifelse(is.na(final$drop_step), .Machine$integer.max, final$drop_step)
  ord <- order(-key, -final$rho1)
  sel <- cand[ord]
  if (!is.null(d)) sel <- sel[seq_len(min(d, length(sel)))]
  
  list(order = sel, selected = cand[final$active], size = length(final$active),
       kappa = final$kappa, threshold = final$threshold, m_reach = final$m_reach,
       truncated = final$truncated, condition = condition,
       init_active = init_active)
}

## -----------------------------------------------------------------------------
## run_method(M, dat, d) -- dispatches to screen_bench()/screen_dpd()/
## screen_tpcc() by M$kind, times it, and combines screening_metrics()
## with runtime and the proxy-whitening errors (err_Sigma,
## err_V), which screen_dpd() computes internally and returns, correctly NA
## for bench/tpcc methods (neither whitens via a proxy).
## -----------------------------------------------------------------------------
run_method <- function(M, dat, d, size_grid = size_grid) {
  p  <- ncol(dat$X) - 1L
  t0 <- proc.time()
  res <- switch(M$kind,
                "bench" = try(screen_bench(dat$y, dat$X, dat$Z, dat$ID,
                                           method = M$method, gamma = M$gamma %||% 0.3),
                              silent = TRUE),
                "dpd"   = try(screen_dpd(dat$y, dat$X, dat$Z, dat$ID, alpha = M$alpha,
                                         proxy = M$proxy, truth = dat$truth,
                                         weighted = M$weighted %||% FALSE),
                              silent = TRUE),
                "dpd_lm"   = try(screen_dpd_lm(dat$y, dat$X, alpha = M$alpha,
                                         truth = dat$truth,
                                         weighted = M$weighted %||% FALSE),
                              silent = TRUE),
                "tpcc"  = try(screen_tpcc(dat$y, dat$X, dat$Z, dat$ID, d = d,
                                          truth = dat$truth), silent = TRUE),
                stop("unknown method kind: ", M$kind))
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  
  failed <- inherits(res, "try-error") || is.null(res$order)
  ord <- if (failed) integer(0) else res$order
  met <- screening_metrics(ord, dat$active, d, p, size_grid = size_grid)
  
  c(met, time = elapsed,
    err_Sigma = if (!failed && !is.null(res$err_Sigma)) res$err_Sigma else NA_real_,
    err_V     = if (!failed && !is.null(res$err_V))     res$err_V     else NA_real_,
    n_nonconv = if (!failed && !is.null(res$n_nonconv)) res$n_nonconv else NA_real_,
    failed = as.numeric(failed))
}


## -----------------------------------------------------------------------------
## screening_metrics(order, active, d, p)
##   order  ranked X-column indices returned by a screening routine
##   active true active screened columns (X-column space, intercept excluded)
##   d      target model size
##
##  TPR      proportion of true signals inside the top d
##  EmpSSP   1 if all signals are inside the top d           (the paper's PIC)
##  MinMS    smallest model size containing every signal = max rank of a signal
##  FP       false positives inside the top d
##  FDP      FP / d
##  rank_*   median and 90th percentile of the ranks of the true signals, and
##           the rank of the worst-ranked signal
##  TPR_at   TPR on a grid of model sizes, for TPR-vs-size curves
## -----------------------------------------------------------------------------
screening_metrics <- function(order, active, d, p,
                              size_grid = seq(10,200, 20)) {
  s <- length(active)
  rk <- match(active, order)                     # NA if never selected
  rk_full <- ifelse(is.na(rk), p, rk)            # censor at p for summaries
  top <- order[seq_len(min(d, length(order)))]
  tp  <- sum(active %in% top)
  out <- c(
    TPR    = tp / s,
    EmpSSP = as.numeric(tp == s),
    MinMS  = max(rk_full),
    FP     = min(d, length(order)) - tp,
    FDP    = (min(d, length(order)) - tp) / max(1, min(d, length(order))),
    rank_med = stats::median(rk_full),
    rank_q90 = as.numeric(stats::quantile(rk_full, 0.9)),
    rank_max = max(rk_full)
  )
  for (g in size_grid)
     out[paste0("TPR_at_", g)] <- sum(active %in% order[seq_len(min(g, length(order)))]) / s
  # out[paste0("TPR_at_", g)] <- sum(active %in% which(beta > n^(-g))) / s

  out
}
