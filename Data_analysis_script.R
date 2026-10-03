## =============================================================================
##  Unified Pipeline for ADNI Real Data Assembly, Screening, and Stability Selection
##  Reproducing the results from the manuscript Ghosh and Thoresen (2026)
## =============================================================================

rm(list = ls())

# ------------------------------------------------------------------------------
# 0. USER CONFIGURATION & GLOBAL PARAMETERS
# ------------------------------------------------------------------------------
target_study   <- "ADNI2"                                  # Target ADNI study cohort
min_visits     <- 3                                        # Minimum visits per subject
cond_vars      <- c("apoe4", "EDUC", "DX", "AGE", "SEX")   # Clinical conditioning covariates
response_var   <- "MMSCORE"                                # Longitudinal response (MMSE score)

# Gene expression preprocessing options
batch_correction <- TRUE                                  # Apply ComBat batch correction
filter_prob      <- TRUE                                  # Apply probe quality filtering
var_quantile     <- 0.00                                  # Lower variance quantile filter cutoff

# File paths
gene_expr_path <- "awData/ADNI_gene_expression_profile.csv"
manifest_path  <- "rawData/gene_probe_manifest.tsv"

# Stability selection settings
run_stability  <- TRUE
B_splits       <- 25                                       # CPSS subsample splits (2B half-runs)
base_seed      <- 20260101L

# Construct baseline output filename prefix
data_file_prefix <- paste0("DerivedData_", target_study, "_", response_var, 
                           "_minVisit", min_visits, "_condVar-", 
                           paste(cond_vars, collapse = "-"))

# Initialize execution logging
log_con <- file(paste0(data_file_prefix, "_pipeline_log.txt"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat(sprintf("=================================================================\n"))
cat(sprintf("ADNI REAL DATA ANALYSIS PIPELINE\n"))
cat(sprintf("Target Cohort: %s | Response: %s | Min Visits: %d\n", target_study, response_var, min_visits))
cat(sprintf("=================================================================\n\n"))

# Load custom functions from local R module folder
if (dir.exists("R")) {
  for (f in list.files("R", pattern = "\\.R$", full.names = TRUE)) source(f)
}


# ------------------------------------------------------------------------------
# 1. DATA ASSEMBLY & PREPROCESSING
# ------------------------------------------------------------------------------
cat("\n--- STEP 1: Assembling Longitudinal Clinical Data ---\n")

suppressPackageStartupMessages(library(ADNIMERGE2))

# Extract repeated MMSE scores from ADQS
qs <- ADNIMERGE2::ADQS
data_long <- qs[qs$PARAMCD == response_var & !is.na(qs$AVAL),
                c("USUBJID", "AVISITN", "VISIT", "AVAL", "AGE", "SEX", "DX", "RACE")]
data_long$RID  <- ADNIMERGE2::convert_usubjid_to_rid(data_long$USUBJID)
data_long$MMSE <- data_long$AVAL
data_long$time <- ADNIMERGE2::convert_number_days(data_long$AVISITN, unit = "year")
data_long      <- data_long[order(data_long$RID, data_long$time), ]

# Compute visit frequencies
n_visits <- table(data_long$RID)
data_long$n_visits  <- as.integer(n_visits[as.character(data_long$RID)])
data_long$visit_num <- ave(data_long$RID, data_long$RID, FUN = seq_along)
rid_visits <- data.frame(RID = as.integer(names(n_visits)), n_visits = as.integer(n_visits))

# Extract APOE4 allele counts from ADSL and GF
adsl <- ADNIMERGE2::ADSL
adsl$RID <- ADNIMERGE2::convert_usubjid_to_rid(adsl$USUBJID)

gf <- ADNIMERGE2::GF
gf <- gf[gf$GFTESTCD == "APOE", c("USUBJID", "GFSTRESC")]
genotype <- gf$GFSTRESC[match(adsl$USUBJID, gf$USUBJID)]
adsl$apoe4 <- ifelse(is.na(genotype), NA_integer_, lengths(regmatches(genotype, gregexpr("4", genotype))))

# Merge clinical covariates into long-format dataset
shared_upper <- setdiff(intersect(toupper(names(data_long)), toupper(names(adsl))), "RID")
dl_names     <- names(data_long)[match(shared_upper, toupper(names(data_long)))]
cl_names     <- names(adsl)[match(shared_upper, toupper(names(adsl)))]

has_match    <- data_long$RID %in\% adsl$RID
clin_matched <- adsl[match(data_long$RID, adsl$RID), cl_names, drop = FALSE]
is_identical <- mapply(function(dl, cl) isTRUE(all.equal(data_long[[dl]][has_match], clin_matched[[cl]][has_match])), dl_names, cl_names)

cols_to_drop <- cl_names[is_identical]
cols_to_add  <- setdiff(names(adsl), c("RID", cols_to_drop))
data_long    <- cbind(data_long, adsl[match(data_long$RID, adsl$RID), cols_to_add, drop = FALSE])

# Filter target study cohort and minimum visit criteria
keep_rid  <- data_long$RID[data_long$ORIGPROT == target_study]
data_long <- data_long[data_long$RID %in% keep_rid, ]

eligible_rids <- rid_visits$RID[rid_visits$n_visits >= min_visits]
data_long     <- data_long[as.character(data_long$RID) %in% eligible_rids, ]
adsl          <- adsl[adsl$RID %in% as.integer(eligible_rids), ]

# Complete cases filter on conditioning variables
rid_na <- unique(data_long$RID[Reduce(`|`, lapply(data_long[cond_vars], is.na))])
data_long <- data_long[!(data_long$RID %in% rid_na), ]

cat(sprintf("Participants retained after clinical filtering: %d (%d longitudinal records)\n", 
            length(unique(data_long$RID)), nrow(data_long)))

cat("\n--- STEP 2: Processing Blood Gene Expression Data ---\n")

# Read gene expression profile CSV
ADNI_METADATA_ROW_LABELS  <- c("Phase", "Visit", "SubjectID", "260/280", "260/230", "RIN", "Affy Plate", "YearofCollection")
ADNI_ANNOTATION_COL_NAMES <- c("ProbeSet", "LocusLink", "Symbol", "Description")
BASELINE_VISIT            <- c("bl", "v03")

raw <- utils::read.csv(gene_expr_path, check.names = FALSE, colClasses = "character", header = FALSE, na.strings = character(0))
is_metadata_row <- trimws(raw[[1]]) %in% ADNI_METADATA_ROW_LABELS
header_row      <- which(trimws(raw[[1]]) == "ProbeSet")[1]
header          <- as.character(raw[header_row, ])
ann_end         <- max(which(trimws(header) %in% ADNI_ANNOTATION_COL_NAMES))
sample_cols     <- seq(ann_end + 1, ncol(raw))

meta_rows   <- raw[is_metadata_row, , drop = FALSE]
sample_meta <- as.data.frame(t(meta_rows[, sample_cols, drop = FALSE]), stringsAsFactors = FALSE)
names(sample_meta) <- trimws(meta_rows[[1]])
sample_meta$RID    <- as.integer(sub("^[0-9]{3}_S_", "", sample_meta$SubjectID))

probe_rows <- raw[(header_row + 1):nrow(raw), , drop = FALSE]
probe_ids  <- probe_rows[[1]]
expr       <- as.matrix(probe_rows[, sample_cols, drop = FALSE])
storage.mode(expr) <- "numeric"
rownames(expr) <- probe_ids
colnames(expr) <- sample_meta$RID

# Gene probe manifest mapping
manifest <- utils::read.delim(manifest_path, colClasses = "character", na.strings = character(0))
manifest$Symbol1 <- vapply(strsplit(manifest$Symbol, " \\|\\| "), `[`, character(1), 1)
manifest$Symbol1[manifest$Symbol == ""] <- NA_character_

# Batch Correction via ComBat
bc_log <- NA
if (batch_correction) {
  cat("Applying sva::ComBat batch correction...\n")
  batch_var <- if ("Affy Plate" %in% names(sample_meta)) "Affy Plate" else NA_character_
  if (!is.na(batch_var) && requireNamespace("sva", quietly = TRUE)) {
    batch <- factor(sample_meta[[batch_var]][match(colnames(expr), sample_meta$RID)])
    keep_batch <- !is.na(batch) & table(batch)[batch] >= 2
    expr_c <- sva::ComBat(dat = expr[, keep_batch, drop = FALSE], batch = droplevels(batch[keep_batch]), par.prior = TRUE)
    bc_log <- list(corrected = TRUE, batch_var = batch_var)
    expr   <- expr_c
  }
}

# Align samples across clinical records and baseline gene expression
common_rid <- intersect(colnames(expr), data_long$RID)
final_rid  <- setdiff(sample_meta$RID[sample_meta$Visit %in% BASELINE_VISIT], rid_visits$RID[rid_visits$n_visits < min_visits])
final_rid  <- intersect(final_rid, common_rid)

data_long <- data_long[as.character(data_long$RID) %in% final_rid, ]
clin      <- adsl[adsl$RID %in% as.integer(final_rid), ]
expr      <- expr[, as.character(final_rid), drop = FALSE]

# Probe Quality & Variance Filtering
n0 <- nrow(expr)
if (filter_prob) {
  # Drop AFFX controls
  expr <- expr[!grepl("^AFFX", rownames(expr)), , drop = FALSE]
  n1   <- nrow(expr)
  # Drop unmapped probes
  expr <- expr[rownames(expr) %in% manifest$ProbeSet[!is.na(manifest$Symbol1)], , drop = FALSE]
  n2   <- nrow(expr)
  # Drop low variance probes
  if (var_quantile > 0) {
    v   <- apply(expr, 1, stats::var, na.rm = TRUE)
    thr <- stats::quantile(v, var_quantile, na.rm = TRUE)
    expr <- expr[v >= thr, , drop = FALSE]
  }
  n3 <- nrow(expr)
  filter_log <- data.frame(step = c("start", "drop_AFFX", "drop_unmapped", "var_quantile"), remaining = c(n0, n1, n2, n3))
}

# Gene-Level Expression Aggregation (select probe with maximum IQR per gene symbol)
sym  <- manifest$Symbol1[match(rownames(expr), manifest$ProbeSet)]
keep <- !is.na(sym)
expr_mapped <- expr[keep, , drop = FALSE]
sym_mapped  <- sym[keep]

iqr       <- apply(expr_mapped, 1, stats::IQR, na.rm = TRUE)
best_idx  <- tapply(seq_along(sym_mapped), sym_mapped, function(idx) idx[which.max(iqr[idx])])
expr_gene <- expr_mapped[unlist(best_idx), , drop = FALSE]
rownames(expr_gene) <- as.character(names(best_idx))

cat(sprintf("Final Preprocessed Dimensions: %d Probes x %d Subjects (%d Genes)\n", 
            nrow(expr), ncol(expr), nrow(expr_gene)))

# Build analysis object and save RDS
ID <- as.integer(factor(data_long$RID))
Cond_data <- model.matrix(~ ., data = data_long[, cond_vars])[, -1]

dat_obj <- list(
  y = data_long$MMSE, ID = ID, Cond_data = Cond_data,
  n = nrow(data_long), m = length(unique(ID)), n_visits = rid_visits,
  expr_probe = expr, expr_gene = expr_gene, gene_meta = sample_meta,
  gene_map = manifest, Clinical_data = clin, data_long = data_long,
  log_info = list(probe_filter = filter_log, batch = bc_log)
)
rds_path <- paste0(data_file_prefix, ".rds")
saveRDS(dat_obj, file = rds_path)
cat(sprintf("Assembly complete. Saved to: %s\n", rds_path))


# ------------------------------------------------------------------------------
# 2. INDEPENDENCE SCREENING (FULL SAMPLE)
# ------------------------------------------------------------------------------
cat("\n--- STEP 3: Executing High-Dimensional Feature Screening ---\n")

# Construct design matrices X and Z
Xp       <- t(dat_obj$expr_probe[, as.character(dat_obj$data_long$RID), drop = FALSE])
X        <- cbind(1, dat_obj$Cond_data, Xp)
colnames(X) <- c("(Intercept)", colnames(dat_obj$Cond_data), rownames(dat_obj$expr_probe))
cond_idx    <- 2:(ncol(dat_obj$Cond_data) + 1)

Z        <- cbind(1, dat_obj$data_long$time)
colnames(Z) <- c("Z0", "Z1")

# Clean missing response/time records if present
keep_rec <- stats::complete.cases(dat_obj$y, Z)
y_clean  <- dat_obj$y[keep_rec]
ID_clean <- dat_obj$ID[keep_rec]
X_clean  <- X[keep_rec, , drop = FALSE]
Z_clean  <- Z[keep_rec, , drop = FALSE]

# Screen size conventions
d_visits  <- target_size(length(y_clean))     # floor(n / log n)
d_subject <- target_size(dat_obj$m)          # floor(m / log m)
dval      <- d_visits

cat(sprintf("Screening target size (d_visits = floor(n / log n)): %d\n", dval))

# Screening methods configuration
methods <- list(
  list(label = "MLE",         kind = "bench",  method = "MLE"),
  list(label = "REML",        kind = "bench",  method = "REML"),
  list(label = "TPCc",        kind = "tpcc"),
  list(label = "cv-P (0.1)",  kind = "dpd",    alpha = 0.1, proxy = "cvP"),
  list(label = "cv-P (0.3)",  kind = "dpd",    alpha = 0.3, proxy = "cvP"),
  list(label = "cv-P (0.5)",  kind = "dpd",    alpha = 0.5, proxy = "cvP"),
  list(label = "I0-P (0.1)",  kind = "dpd",    alpha = 0.1, proxy = "I0"),
  list(label = "I0-P (0.3)",  kind = "dpd",    alpha = 0.3, proxy = "I0"),
  list(label = "I0-P (0.5)",  kind = "dpd",    alpha = 0.5, proxy = "I0"),
  list(label = "no-P (0.1)",  kind = "dpd_lm", alpha = 0.1, proxy = "none"),
  list(label = "no-P (0.3)",  kind = "dpd_lm", alpha = 0.3, proxy = "none"),
  list(label = "no-P (0.5)",  kind = "dpd_lm", alpha = 0.5, proxy = "none")
)

screening_results <- list()

for (M in methods) {
  cat(sprintf("Running Screening: %-12s (d = %d, cond_cols = %d)... ", M$label, dval, length(cond_idx)))
  t0 <- proc.time()
  
  res <- switch(M$kind,
                "bench"  = try(screen_bench(y_clean, X_clean, Z_clean, ID_clean, Cond = cond_idx, method = M$method, gamma = M$gamma %||% 0.3), silent = TRUE),
                "dpd"    = try(screen_dpd(y = y_clean, X = X_clean, Z = Z_clean, ID = ID_clean, alpha = M$alpha, Cond = cond_idx, proxy = M$proxy, weighted = M$weighted %||% FALSE), silent = TRUE),
                "dpd_lm" = try(screen_dpd_lm(y_clean, X_clean, alpha = M$alpha, Cond = cond_idx, weighted = M$weighted %||% FALSE), silent = TRUE),
                "tpcc"   = try(screen_tpcc(y_clean, X_clean, Z_clean, ID_clean, d = dval, Cond = cond_idx), silent = TRUE),
                stop("Unknown method kind: ", M$kind)
  )
  
  elapsed <- as.numeric((proc.time() - t0)["elapsed"])
  
  if (inherits(res, "try-error") || is.null(res$order)) {
    cat("FAILED.\n")
    next
  }
  
  top <- res$order[seq_len(min(dval, length(res$order)))]
  gene_cols <- top[top > max(cond_idx)]
  
  screening_results[[paste0(M$label, "_d", dval)]] <- list(
    method  = M$label,
    d       = dval,
    runtime = elapsed,
    probes  = colnames(X_clean)[gene_cols],
    genes   = dat_obj$gene_map$Symbol1[match(colnames(X_clean)[gene_cols], dat_obj$gene_map$ProbeSet)]
  )
  cat(sprintf("Done (%.2f s).\n", elapsed))
}

screening_file <- paste0(data_file_prefix, "_screeningResults.rds")
saveRDS(screening_results, screening_file)

if (file.exists("export_wide.R")) {
  source("export_wide.R")
  export_all_wide(screening_results, outdir = "summary_wide")
}

if (file.exists("summary_results.R")) {
  source("summary_results.R")
  summarize_screening(screening_results, outdir = "summary_screening")
}


# ------------------------------------------------------------------------------
# 3. COMPLEMENTARY PAIR STABILITY SELECTION (CPSS)
# ------------------------------------------------------------------------------
if (run_stability) {
  cat("\n--- STEP 4: Executing Stability Selection (CPSS) ---\n")
  
  if (file.exists("stability_selection.R") && file.exists("build_gene_stability_table.R")) {
    source("stability_selection.R")
    source("build_gene_stability_table.R")
    
    cpss_out <- vector("list", length(methods))
    names(cpss_out) <- vapply(methods, `[[`, character(1), "label")
    
    for (i in seq_along(methods)) {
      M <- methods[[i]]
      cat(sprintf("CPSS: %-12s (d = %d, B = %d -> %d half-runs)... ", M$label, dval, B_splits, 2L * B_splits))
      t0 <- proc.time()
      
      cpss_res <- try(cpss_stability(M, y_clean, X_clean, Z_clean, ID_clean, cond_idx, dval, B = B_splits, seed = base_seed + i), silent = TRUE)
      
      if (inherits(cpss_res, "try-error")) {
        cat("FAILED.\n")
        cpss_out[[i]] <- NULL
      } else {
        elapsed <- as.numeric((proc.time() - t0)["elapsed"])
        cat(sprintf("Done (%.1f s, %d/%d valid half-runs).\n", elapsed, cpss_res$n_success, cpss_res$n_runs))
        cpss_out[[i]] <- cpss_res
      }
    }
    
    cpss_clean <- cpss_out[!vapply(cpss_out, is.null, logical(1))]
    saveRDS(cpss_clean, file = paste0(data_file_prefix, "_cpss_results.rds"))
    
    if (exists("summarize_stability")) {
      summarize_stability(cpss_clean, X_clean, dat_obj$gene_map,
                          export_probes = screening_file,
                          export_genes  = screening_file)
    }
    
    if (exists("build_gene_stability_table")) {
      build_gene_stability_table(screening_file, paste0(data_file_prefix, "_cpss_results.rds"), dval)
    }
  } else {
    cat("Stability selection modules (stability_selection.R / build_gene_stability_table.R) not found. Skipping.\n")
  }
}

cat("\n=================================================================\n")
cat("Pipeline execution completed successfully.\n")
cat("=================================================================\n")

sink(type = "output")
sink(type = "message")