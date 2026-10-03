# ==============================================================================
# Simulation Study for Ghosh and Thoresen (2026): Reproducibility of Manuscript Results
# ==============================================================================
#
# Description:
#   This script executes Monte Carlo simulation studies to evaluate screening 
#   performance under various simulation specification and contamination schemes,
#   as presented in Ghosh and Basu (2026) paper 
#
# Inputs:
#   - User configuration settings defined at the top of this script (Section 1).
#   - Helper functions in the 'R/' directory (automatically sourced).
#
# Outputs:
#   - Raw CSV result files saved to 'Results (Raw)/'
#   - Summary statistics CSV tables saved to working directory
#   - Diagnostic and publication-ready plots saved to 'Figures/'
# ==============================================================================

# Clear workspace
rm(list = ls())

# ------------------------------------------------------------------------------
# Section 1 (USER CONFIGURATION)
# ------------------------------------------------------------------------------

# Simulation dimensions and replication counts
# Execution Mode: specify p (number of FE covariates) and ITS (number of MC replications)
p <- 1000;  ITS <- 100
q     <- 4

# Design Matrix Strategy:
#   FIXED_DESIGN = TRUE  -> Generate X and Z once and reuse across replications
#   FIXED_DESIGN = FALSE -> Regenerate X and Z in every replication
FIXED_DESIGN <- TRUE

# Parallel Processing:   Number of CPU cores to leave unused for system stability
LeaveCore <- 1

# Random Seed for Reproducibility:
BASE_SEED <- 20260101L

# Random-effects structure ("R1" or "R2") as per the manuscript:
re <- "R1"

# Variance of random effects (sigma2_b) and errors (sigma)
sigma2_b <- 1
sigma <- 1

# Random-effects correlation parameter:
rho_b <- 0.3

# True Signal Sparsity Pattern as per the manuscript:
#   "S1" -> Active predictors at positions 1:5
#   "S2" -> Active predictors at positions c(5, 10, 50, 100, 200)
signal <- "S2"

# Signal decay parameter (signal = n^(-kappa)):
#   kappa = 0   : Signal magnitude constant with sample size
#   kappa = 0.2 : Signal magnitude decreases with sample size
kappa <- 0.2


# Covariance structure of X:
#   "I"   -> Independent predictors
#   "CS"  -> Compound Symmetry
#   "AR1" -> First-order Autoregressive (Toeplitz)
sigx <- "CS"

# Fixed-effects predictor correlation (for CS/AR1 structures):
rho_x <- 0.5


# Cluster Design:
m  <- 10                              # Number of clusters
ni <- rep(c(10, 15, 20, 25, 30), 2)   # Observations per cluster (vector or scalar)


# Error and Random-Effect Distributions:
#   "gauss" -> Normal errors + Normal random effects
#   "heavy" -> Heavy-tailed t(df=5) errors + Normal random effects
#   "skew"  -> Skew-normal errors + Normal random effects
#   "RE_t5" -> Normal errors + Multivariate-t(df=5) random effects
dist <- "gauss"

# Grid of retained model sizes evaluated during variable selection:
sizes <- seq(10, 200, 20)


# ------------------------------------------------------------------------------
# 2. PACKAGE DEPENDENCIES & ENVIRONMENT SETUP
# ------------------------------------------------------------------------------

# Required and optional packages
.req_pkgs <- c("MASS", "utils", "mvtnorm", "lme4", "glmnet", 
               "Matrix", "foreach", "doParallel", "fsimR")
.opt_pkgs <- c("doRNG", "sn", "robustlmm", "RhpcBLASctl")

# Install and load missing required packages
for (p in .req_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    message("Installing missing package: ", p)
    install.packages(p, repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(p, character.only = TRUE))
}

# Load optional packages if available
for (p in .opt_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    suppressPackageStartupMessages(library(p, character.only = TRUE))
  }
}

# Define null-coalescing operator if not already loaded
if (!exists("%||%")) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

# Source all custom R functions from the 'R' directory
if (dir.exists("R")) {
  r_files <- list.files("R", pattern = "\\.R$", full.names = TRUE)
  sapply(r_files, source)
} else {
  stop("Directory 'R/' containing helper functions was not found.")
}

# ------------------------------------------------------------------------------
# 3. PARALLEL COMPUTING SETUP
# ------------------------------------------------------------------------------

# Limit multi-threading in BLAS/LAPACK backends to avoid thread contention
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

# Initialize parallel cluster
nprocs  <- max(1L, parallel::detectCores() - LeaveCore)
para_cl <- parallel::makeCluster(nprocs)
doParallel::registerDoParallel(para_cl)

# Ensure cluster cleanup on script exit or interruption
on.exit({
  parallel::stopCluster(para_cl)
  message("Parallel cluster successfully terminated.")
}, add = TRUE)

# ------------------------------------------------------------------------------
# 4. DERIVED SIMULATION PARAMETERS & METHOD SPECIFICATIONS
# ------------------------------------------------------------------------------

# Cluster structure calculations
cl <- list(m = m, ni = ni)
ni <- if (length(cl$ni) == 1) rep(cl$ni, cl$m) else cl$ni
n  <- sum(ni)


# Signal strength decaying factor
B_val <- n^(-kappa)

# Distribution specification map
dstr <- switch(dist,
               "gauss" = list(err_dist = "norm", re_dist = "mvnorm"),
               "heavy" = list(err_dist = "t",    re_dist = "mvnorm", df = 5),
               "skew"  = list(err_dist = "sn",   re_dist = "mvnorm", alpha = 5),
               "RE_t5" = list(err_dist = "norm", re_dist = "mvt",    df = 5),
               stop("Unknown distribution choice: ", dist)
)

# Consolidated configuration list
cfg <- c(list(
  m = cl$m, ni = ni, p = p, q = q, signal = signal, re = re,
  sigx = sigx, rho_x = rho_x, rho_b = rho_b, sigma = sigma,
  B = B_val, Bspec = Bspec, clus = paste0("m", cl$m, "_ni", paste(unique(ni), collapse = "-")),
  dist = dist, n = n, d = target_size(n)
), dstr)

# Evaluation methods grid
meths <- list(
  # Benchmarks
  list(label = "MLE",          kind = "bench", method = "MLE"),
  list(label = "REML",         kind = "bench", method = "REML"),
  list(label = "TPCc",         kind = "tpcc"),
  
  # Oracle-proxy Density Power Divergence (DPD)
  list(label = "orcl-P (0.1)", kind = "dpd",   alpha = 0.1, proxy = "oracle"),
  list(label = "orcl-P (0.3)", kind = "dpd",   alpha = 0.3, proxy = "oracle"),
  list(label = "orcl-P (0.5)", kind = "dpd",   alpha = 0.5, proxy = "oracle"),
  
  # Cross-validated proxy DPD
  list(label = "cv-P (0.1)",   kind = "dpd",   alpha = 0.1, proxy = "cvP"),
  list(label = "cv-P (0.3)",   kind = "dpd",   alpha = 0.3, proxy = "cvP"),
  list(label = "cv-P (0.5)",   kind = "dpd",   alpha = 0.5, proxy = "cvP"),
  
  # Identity proxy DPD
  list(label = "I0-P (0.1)",   kind = "dpd",   alpha = 0.1, proxy = "I0"),
  list(label = "I0-P (0.3)",   kind = "dpd",   alpha = 0.3, proxy = "I0"),
  list(label = "I0-P (0.5)",   kind = "dpd",   alpha = 0.5, proxy = "I0"),
  
  # No-proxy DPD
  list(label = "no-P (0.1)",   kind = "dpd",   alpha = 0.1, proxy = "none"),
  list(label = "no-P (0.3)",   kind = "dpd",   alpha = 0.3, proxy = "none"),
  list(label = "no-P (0.5)",   kind = "dpd",   alpha = 0.5, proxy = "none")
)

# Contamination schemes definition
CONTAM <- data.frame(
  scheme = c("C0", "C1", "C1", "C2", "C2", "C3", "C3", "C4", "C4", "C5"),
  prop   = c(0,    0.05, 0.10, 0.05, 0.10, 0.05, 0.10, 0.05, 0.10, 0.20),
  stringsAsFactors = FALSE
)

CONTAM_TARGET <- c(C2 = "active", C3 = "inactive")

# ------------------------------------------------------------------------------
# 5. LOGGING & DATA GENERATION
# ------------------------------------------------------------------------------

# Build output tag and initialize log file
tag <- sprintf("%s_%s_%s(%s)_RE(%s_%s)_%s_%s_p%d_Bk%s",
               cfg$signal, cfg$re, cfg$sigx, rho_x,
               sigma2_b, rho_b, cfg$clus, cfg$dist, cfg$p, 10 * kappa)

log_file <- sprintf("%s_%s.log", tag, format(Sys.time(), "%Y%m%d-%H%M%S"))
cat(sprintf("# Log opened %s\n", Sys.time()), file = log_file)

logf <- function(fmt, ...) {
  msg <- sprintf(fmt, ...)
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), msg),
      file = log_file, append = TRUE)
  invisible(msg)
}

logf("Simulation initialized: n=%d, d=%d, methods=%d", cfg$n, cfg$d, length(meths))

# Setup true parameter values and covariance matrices
pos <- if (cfg$signal == "S1") 1:5 else c(5, 10, 50, 100, 200)
pos <- pos[pos <= p]

beta <- numeric(p)
beta[pos] <- cfg$B
active <- pos
is_R1 <- identical(cfg$re, "R1")

dm <- target_size(m)
sizes <- unique(sort(c(dm - 1, dm, dm + 3, sizes)))
sizes <- sizes[sizes <= cfg$p]

# Covariance matrix for X
SigmaX <- construct_covMat(
  p,
  type = switch(cfg$sigx,
                "I"   = "diag",
                "CS"  = "CS",
                "AR1" = "toeplitz",
                stop("Unknown sigx type: ", cfg$sigx)),
  params = list(rho = cfg$rho_x, variances = 1)
)

# Covariance matrix for Random Effects
Psi <- construct_covMat(q, type = "CS", corr = FALSE, params = list(rho = rho_b, variances = sigma2_b))
diag(Psi) <- sigma2_b

# Distribution parameters bundle
dist_settings <- list(
  X_distr = list(
    distr_name   = "mvnorm",
    distr_params = list(sigma = SigmaX)
  ),
  RE_distr = list(
    distr_name   = cfg$re_dist %||% "mvnorm",
    distr_params = distr_spec(sigma = Psi, df = cfg$df %||% 5, alpha = cfg$alpha %||% 5)
  ),
  error_distr = list(
    distr_name   = cfg$err_dist %||% "norm",
    distr_params = list(mean = 0, sigma = cfg$sigma, df = cfg$df  %||% 5, alpha = cfg$alpha %||% 5)
  ),
  Z_distr = if (is_R1) NULL else list(distr_name = "mvnorm", distr_params = list(dim = q - 1))
)

sim_settings <- sim_spec(
  n_subj = cfg$m, n_obs = ni, p = p, q = q,
  beta_coeff = beta, is.ZInX = is_R1, include.Xintercept = FALSE
)

# Generate baseline uncontaminated data across replications
n_rep <- if (FIXED_DESIGN) c(1, ITS) else c(ITS, 1)
all_data <- simulate_LMMdata(n_rep = n_rep, sim_settings = sim_settings,
                             distr_settings = dist_settings, seed = BASE_SEED)

truth <- list(Psi = Psi, sigma2 = cfg$sigma^2)
export_vars <- ls(envir = .GlobalEnv)

# ------------------------------------------------------------------------------
# 6. MAIN SIMULATION LOOP (CONTAMINATION SCHEMES)
# ------------------------------------------------------------------------------

all_res <- list()
if (!dir.exists("Results (Raw)")) dir.create("Results (Raw)", recursive = TRUE)

for (ci in seq_len(nrow(CONTAM))) {
  scheme <- CONTAM$scheme[ci]
  prop   <- CONTAM$prop[ci]
  target <- if (scheme %in% names(CONTAM_TARGET)) CONTAM_TARGET[[scheme]] else "inactive"
  
  out_file <- file.path("Results (Raw)", sprintf("%s_%s_prop%02d.csv", tag, scheme, round(100 * prop)))
  
  if (file.exists(out_file)) {
    logf("Skip existing output: %s", basename(out_file))
    all_res[[ci]] <- utils::read.csv(out_file, stringsAsFactors = FALSE)
    next
  }
  
  # Apply contamination scheme
  sim_data <- all_data
  
  if (scheme == "C1") {
    sim_data <- add_contamination(
      sim_data, cont_pos = list(name = "y"),
      cont_settings = list(
        cont_mode = "casewise", cont_prop = prop, cont_type = "additive",
        cont_distr = list(distr_name = "norm", distr_params = list(mean = 20, sd = 1))
      )
    )
  } else if (scheme == "C2") {
    col <- pick_col(target, active, p)
    sim_data <- add_contamination(
      sim_data, cont_pos = list(name = "X", col = col),
      cont_settings = list(
        cont_mode = "casewise", cont_prop = prop, cont_type = "additive",
        cont_distr = list(distr_name = "norm", distr_params = list(mean = -2, sd = 1))
      )
    )
  } else if (scheme == "C3") {
    col <- pick_col(target, active, p)
    sim_data <- add_contamination(
      sim_data, cont_pos = list(name = "X", col = col),
      cont_settings = list(
        cont_mode = "casewise", cont_prop = prop, cont_type = "leverage",
        leverage_factor = 5, outlier_factor = 5
      )
    )
  } else if (scheme == "C4") {
    stopifnot(ncol(sim_data$Z) >= 2)
    sim_data <- add_contamination(
      sim_data, cont_pos = list(name = "Z", col = 2),
      cont_settings = list(
        cont_mode = "casewise", cont_prop = prop, cont_type = "leverage",
        leverage_factor = 5, outlier_factor = 5
      )
    )
  } else if (scheme == "C5") {
    m_clusters <- length(unique(sim_data$ID))
    nc <- max(1L, round(prop * m_clusters))
    sel_clusters <- sample(unique(sim_data$ID), nc)
    rows <- which(sim_data$ID %in% sel_clusters)
    
    sub <- list(y = sim_data$y[rows], X = sim_data$X[rows, , drop = FALSE], ID = sim_data$ID[rows])
    sub <- add_contamination(
      sub, cont_pos = list(name = "y"),
      cont_settings = list(
        cont_mode = "casewise", cont_prop = 1, cont_type = "additive",
        cont_distr = list(distr_name = "norm", distr_params = list(mean = 20, sd = 1))
      )
    )
    col <- pick_col("inactive", active, p)
    sub <- add_contamination(
      sub, cont_pos = list(name = "X", col = col),
      cont_settings = list(
        cont_mode = "casewise", cont_prop = 1, cont_type = "replace",
        cont_distr = list(distr_name = "norm", distr_params = list(mean = 0, sd = 5))
      )
    )
    sim_data$y[rows]      <- sub$y
    sim_data$X[rows, col] <- sub$X[, col]
  }
  
  # Execute Monte Carlo replications in parallel
  t_start <- Sys.time()
  res <- foreach(
    r = seq_len(ITS), .combine = "rbind",
    .packages = c("lme4", "mvtnorm", "glmnet", "MASS", "fsimR"),
    .export = export_vars
  ) %dopar% {
    set.seed(BASE_SEED + r)
    
    if (FIXED_DESIGN) {
      dat <- list(y = as.numeric(sim_data$y[, r]), X = cbind(1, sim_data$X),
                  Z = sim_data$Z, ID = as.integer(sim_data$ID),
                  active = active + 1, truth = truth)
    } else {
      dat <- list(y = as.numeric(sim_data$y[[r]]), X = cbind(1, sim_data$X[[r]]),
                  Z = sim_data$Z[[r]], ID = as.integer(sim_data$ID),
                  active = active + 1, truth = truth)
    }
    dat$Z <- standardize_Z(dat$Z)
    
    do.call(rbind, lapply(meths, function(M) {
      v <- run_method(M, dat, cfg$d, size_grid = sizes)
      data.frame(rep = r, method = M$label, as.list(v),
                 cont_col = dat$meta$col %||% NA,
                 cont_col_active = dat$meta$col_active %||% NA,
                 stringsAsFactors = FALSE)
    }))
  }
  
  # Append simulation metadata and export results
  res <- cbind(res, signal = cfg$signal, re = cfg$re, sigx = cfg$sigx,
               rho_x = cfg$rho_x, rho_b = cfg$rho_b, clus = cfg$clus,
               dist = cfg$dist, p = cfg$p, n = cfg$n, d = cfg$d,
               Bspec = cfg$Bspec, B = cfg$B, scheme = scheme, prop = prop)
  
  write.csv(res, out_file, row.names = FALSE)
  all_res[[ci]] <- res
  logf("Completed scheme %s (%s) in %.1fs", scheme, basename(out_file),
       as.numeric(Sys.time() - t_start, units = "secs"))
}

# ------------------------------------------------------------------------------
# 7. SUMMARY TABLES AND DIAGNOSTIC PLOTS
# ------------------------------------------------------------------------------

all_res <- do.call(rbind, all_res)
all_res$cell <- paste(all_res$scheme, sprintf("\%.2f", all_res$prop), sep = "|")
cell_f <- factor(all_res$cell, levels = unique(all_res$cell))

method_levels <- unique(all_res$method)
method_levels <- method_levels[!is.na(method_levels)]

# 1. Main Summary Table
summary_tab <- do.call(rbind, lapply(split(all_res, cell_f), function(x) {
  s <- summarise_runs(x, by = "method")
  cbind(scheme = x$scheme[1], prop = x$prop[1], method = s$key, s[, -1])
}))

summary_tab$paper_cell <- sprintf("%.2f (%.2f)", summary_tab$TPR_med, summary_tab$EmpSSP)
summary_file <- paste0(tag, "_summary.csv")
write.csv(summary_tab, summary_file, row.names = FALSE)
logf("Summary table exported: %s", summary_file)

# 2. Proxy Error Summary Table
proxy_tab <- do.call(rbind, lapply(split(all_res, cell_f), function(x) {
  do.call(rbind, lapply(split(x, x$method), function(xx)
    data.frame(scheme = xx$scheme[1], prop = xx$prop[1], method = xx$method[1],
               errSigma_med = stats::median(xx$err_Sigma, na.rm = TRUE),
               errV_med     = stats::median(xx$err_V,     na.rm = TRUE),
               EmpSSP       = mean(xx$EmpSSP, na.rm = TRUE),
               stringsAsFactors = FALSE)))
}))

proxy_file <- file.path("Results (Raw)", paste0(tag, "_proxy_errors.csv"))
write.csv(proxy_tab, proxy_file, row.names = FALSE)
logf("Proxy errors table exported: %s", proxy_file)

# 3. Create Figures Directory and Plot Graphics
if (!dir.exists("Figures")) dir.create("Figures", recursive = TRUE)

# A. Minimum Model Size Boxplots
pdf_minms <- file.path("Figures", paste0(tag, "_MinMS_boxplots.pdf"))
grDevices::pdf(pdf_minms, width = 12, height = 6)
for (k in unique(all_res$cell)) {
  x  <- all_res[all_res$cell == k, , drop = FALSE]
  mm <- tapply(x$MinMS, list(x$rep, x$method), identity)
  method_levels_k <- method_levels[method_levels %in% colnames(mm)]
  mm <- mm[, method_levels_k, drop = FALSE]
  
  boxplot_enhanced(as.data.frame(mm), xlab = "", title = paste0("Contamination Cell: ", k),
                   ylab = "Minimum Model Size", cex.axis = 1.2, cex.lab = 1.5, font.axis = 2)
  graphics::title(main = paste(tag, k, sep = " | "), line = 3)
}
grDevices::dev.off()
logf("MinMS boxplots saved: %s", pdf_minms)

# B. Runtime Boxplots
pdf_runtime <- file.path("Figures", paste0(tag, "_runtime_boxplots.pdf"))
grDevices::pdf(pdf_runtime, width = 12, height = 6)
for (k in unique(all_res$cell)) {
  x  <- all_res[all_res$cell == k, , drop = FALSE]
  rt <- tapply(x$time, list(x$rep, x$method), identity)
  method_levels_k <- method_levels[method_levels %in% colnames(rt)]
  rt <- rt[, method_levels_k, drop = FALSE]
  
  boxplot_enhanced(as.data.frame(rt), xlab = "", ylab = "Runtime (sec)",
                   cex.axis = 1.2, cex.lab = 1.5, font.axis = 2)
  graphics::title(main = paste(tag, k, sep = " | "), line = 3)
}
grDevices::dev.off()
logf("Runtime boxplots saved: %s", pdf_runtime)

# C. Plot TPR Curves Across Retained Model Sizes (d)
pdf_tpr <- file.path("Figures", paste0(tag, "_TPRcurve.pdf"))
cols <- grDevices::hcl.colors(length(method_levels), palette = "Dark 3")
ltys <- rep(c(1, 2, 3, 4, 5, 6), length.out = length(method_levels))

grDevices::pdf(pdf_tpr, width = 7.5, height = 5.8, useDingbats = FALSE)
old_par <- graphics::par(no.readonly = TRUE)

graphics::par(
  mar = c(4.5, 4.8, 1.5, 1.2), mgp = c(2.7, 0.8, 0),
  tcl = -0.25, las = 1, cex.axis = 1.05, cex.lab = 1.15,
  cex.main = 1.1, font.lab = 1, xaxs = "i", yaxs = "i"
)

for (k in unique(all_res$cell)) {
  df <- all_res[all_res$cell == k, , drop = FALSE]
  
  graphics::plot(
    range(sizes), c(0, 1), type = "n", log = "x",
    xlim = range(sizes), ylim = c(0, 1.1),
    xlab = expression("Retained model size " * italic(d)),
    ylab = "True Positive Rate (TPR)", xaxt = "n"
  )
  
  graphics::axis(1, at = sizes, labels = sizes)
  graphics::abline(h = seq(0, 1, by = 0.1), 
                   col = grDevices::adjustcolor("grey70", alpha.f = 0.35), lty = 1, lwd = 0.7)
  graphics::box()
  graphics::abline(v = c(target_size(n), target_size(cl$m)), col = c("grey35", "grey35"), lty = c(2, 3), lwd = 1.5)
  
  for (i in seq_along(method_levels)) {
    method_i <- method_levels[i]
    dfi <- df[df$method == method_i, , drop = FALSE]
    if (nrow(dfi) == 0) next
    
    tpr_cols <- paste0("TPR_at_", sizes)
    tpr_cols <- tpr_cols[tpr_cols %in% names(dfi)]
    if (length(tpr_cols) == 0) next
    
    tpr_mean <- colMeans(dfi[, tpr_cols, drop = FALSE], na.rm = TRUE)
    sizes_i  <- sizes[paste0("TPR_at_", sizes) %in% tpr_cols]
    
    graphics::lines(sizes_i, tpr_mean, type = "o", pch = 16, cex = 0.65, lwd = 2.2, lty = ltys[i], col = cols[i])
  }
  
  present_methods <- method_levels[method_levels %in% unique(df$method)]
  present_i <- match(present_methods, method_levels)
  
  graphics::legend(
    "bottomright", legend = present_methods, col = cols[present_i],
    lty = ltys[present_i], lwd = 2.2, pch = 16, pt.cex = 0.65, bty = "n",
    bg = grDevices::adjustcolor("white", alpha.f = 0.85), inset = 0.02, cex = 0.9
  )
  
  graphics::mtext(text = paste(tag, k, sep = " | "), side = 3, line = 0.2, cex = 0.95, font = 2)
}

graphics::par(old_par)
grDevices::dev.off()
logf("TPR curves PDF generated: %s", pdf_tpr)

logf("Simulation study execution completed successfully.")