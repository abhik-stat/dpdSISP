# dpdSISP: Robust Sure Screening for Ultrahigh-Dimensional Linear Mixed Models using DPD-SISP

**DPD-SISP** is a minimum density power divergence (DPD) estimator based marginal screening procedure for
fixed-effects in linear mixed models. It first *whitens* the data using a proxy for the unknown
random-effects covariance structure, then ranks candidate covariates by a robust marginal
association statistic, designed to remain reliable when a fraction of observations are
contaminated, and to scale to tens of thousands of candidates.

This repository contains R Codes accompanying:

> Ghosh, A. and Thoresen, M. **Robust and Scalable Sure Screening of Fixed Effects in
> Ultrahigh-dimensional Linear Mixed Models,** which is presently under review in *Statistics and Computing*.

<!-- badges: start -->

<!-- badges: end -->

<br>


---

## General Function Library

Every screening procedures used in our study (Ghosh and Thoresen, 2026) is defined 
as a generic functions, along with necessary background and heper function, in `R/`. 

<br>

### Requirements

`R` (>= 4.1 recommended).

```r
install.packages(c("MASS", "mvtnorm", "lme4", "glmnet", "Matrix"))
```

---
---


### Repository structure

```
R/
  lmm-estimates.R         Benchmark LMM fitters 
  mdpde.R                 Minimum DPD estimation 
  proxy.R                 Random-effects proxies and whitening 
  screening.R             Screening procedures and metrics 
  utilities.R             Shared helpers 
```

### `lmm-estimates.R` -- benchmark LMM fitters

| Function | What it does |
|---|---|
| `LMM_MLE()` | Maximum likelihood estimation under LMM, using `lme4::lmer()` |
| `LMM_REML()` | Restricted maximum likelihood estimation under LMM, using `lme4::lmer()` |
| `LMM_HGD()` | Minimum hierarchical gamma-divergence estimator (Sugasawa et al., 2025)), fit by an MM algorithm|

### `mdpde.R` -- minimum DPD estimation

| Function | What it does |
|---|---|
| `lmdpd()` | Minimum DPD estimator for a linear regression (Ghosh and Basu, 2013), using IRLS or the L-BFGS-B as fallback |
| `mdpde_lmm_icept()` | Cluster-wise minimum DPD estimation of `(sigma^2, Psi)`, following Saraceno et al. (2024), in the intercept-only LMM (used by the rob-P proxy)|


### `proxy.R` -- random-effects proxies and whitening

| Function | What it does |
|---|---|
| `make_proxy()` | Builds the q x q random-effects proxy matrix. Available options are oracle Proxy (truth), cv-P, I0-P, rob-P, none (no proxy), a deliberately misspecified scaling, or the fixed log(n) choice |
| `whiten_by_cluster()` | Applies the block-diagonal whitening transform, cluster by cluster |
| `I0_proxy()` | The intercept-only-REML proxy (I0-P) |
| `cvP_proxy()` | The cross-validated scalar proxy (cv-P), where folds are taken over clusters |
| `standardize_Z()` | Rescales each non-intercept column of `Z` to unit root-mean-square |
| `inv_sqrt_sym()`, `DZ_scale()` | Helpers used by the other functions |
| `proxy_diagnostics()` | A proxy's realized spectral accuracy against the true covariance  (it needs the truth) |


### `screening.R` -- screening procedures and metrics

| Function | What it does |
|---|---|
| `screen_dpd()` | DPD-SISP / DPD-CSISP, the paper's proposed robust screening procedures, based on any proxy from `make_proxy()` and any robustness parameter `alpha` |
| `screen_dpd_lm()` | The same marginal DPD statistic computed directly on the (conditioned) response, with no random-effects term and no whitening at all (equivalent to screen_dpd() with proxy as `none`|
| `screen_bench()` | Benchmarks scrreening, via `lmm-estimates.R`'s `LMM_MLE`/`LMM_REML`/`LMM_HGD` |
| `screen_isis()` | DPD-ISISP, the iterative extension, for settings where one-step marginal screening is insufficient |
| `screen_tpcc()` | TPCc, the conditional thresholded partial correlation (Alabiso & Shang, 2023), applying the TPC algorithm of Li, Liu & Lou (2017) conditional on the random effects (see Notes below) |
| `run_tpc_steps()` | The Step-1-3 engine `screen_tpcc()` calls twice: once unconditioned, once on the BLUP-adjusted response |
| `tpc_kurtosis()`, `tpc_threshold()`, `pcor_given()` | TPC's own excess-kurtosis estimate, rejection threshold, and partial-correlation helper |
| `run_method()` | Dispatches to the right screening function by a method's `kind`, times it, and attaches `screening_metrics()` |
| `screening_metrics()` | TPR, EmpSSP, MinMS, FP, FDP, signal ranks, and TPR at a grid of sizes, from a ranked selection against the true active set |

### `utilities.R` -- shared helpers

| Function | What it does |
|---|---|
| `target_size()` | The target model size, `floor(n / log(n))`, using  natural logarithm |
| `split_clusters()` | Splits cluster IDs into two halves, for an independent proxy-estimation/screening split |
| `pick_col()` | Resolves a contamination target (`"active"`/`"inactive"`/a specific column) to an actual column index |
| `summarise_runs()` | Aggregates replications into the paper's table format, reporting defiiferent summary measures including median TPR,  EmpSSP, etc. |
| `boxplot_enhanced()` | Personalized enhanced boxplot/violin plotting function, used to generate some figures of the main paper |
| `%||%` | Null-coalescing operator (`a %||% b` is `a` unless `NULL`, else `b`)|

### Notes on the implementation

- **Whitening is cluster-wise (block-diagonal), not a full n x n transform.** The working
  covariance under the paper's model is block-diagonal across clusters; `whiten_by_cluster()`
  whitens cluster by cluster, both for statistical correctness (a full-matrix transform mixes
  clusters) and for computational cost.
- **TPCc's citation is two papers**: the underlying TPC algorithm is Li, Liu and Lou
  (2017, *Statistica Sinica* 27, 983-996); applying it conditional on the random effects, as
  `screen_tpcc()` does, is Alabiso and Shang (2023, *Communications in Statistics -- Theory and
  Methods* 52, 6355-6380). `condition = "blup"` (the default) subtracts the BLUPs of an
  intercept-only LMM before applying TPC, matching Alabiso & Shang's Section 3.1; 
  `"whiten"` and `"none"` are our own alternatives, kept for comparison, not from either paper.
- **The `oracle` proxy is a diagnostic upper benchmark**, using the true random-effects
  covariance; it needs a ground truth real data does not have. It is used in simulation, 
  only to examine the cost of proxy estimation.


---


## Reproducing Simulation Study from the Paper

This section provides instructions for reproducing the simulation results presented 
in the manuscript Ghosh and Thoresen (2026) using the master script `simulation_script.R`.

The main script `simulation_script.R` executes Monte Carlo experiments comparing 
various screening procedures under high-dimensional LMM across various 
contamination schemes, (C0)--(C5) as specified in the manuscript, 
for a given simulation specification.

The script automatically:

1. Loads required core packages and helper routines.
2. Configures single-core operations for BLAS/LAPACK to prevent thread oversubscription.
3. Initializes a parallel backend (`doParallel` / `foreach`) reserving standard system CPU cores.
4. Generates baseline LMM datasets using `fsimR` package.
5. Iterates through contamination scenarios, running specified benchmark and DPD methods.
6. Exports raw CSV results, summary tables, and publication-ready PDF plots.

---

### Configuration Parameters

All user-configurable parameters are defined at the beginning of 
`simulation_script.R` under **Section 1 (USER CONFIGURATION)**:

| Configuration Parameter | Default (Specific Options) | Description |
| --- | --- | --- |
| `p` | `1000` | Number of fixed-effect covariates ($p$). |
| `ITS` | `100` | Number of Monte Carlo replications. |
| `q` | `4` | Dimension of the random-effects vector ($q$). |
| `FIXED_DESIGN` | `TRUE` (`TRUE` / `FALSE`) | Design matrix regime. `TRUE` fixes design matrices $X$ and $Z$ across Monte Carlo replications; `FALSE` regenerates design matrices in every replication. |
| `LeaveCore` | `1` (`0`, `1`) | `1` means a CPU core is left unused for system stability during parallel processing. |
| `BASE_SEED` | `20260101L` | Global seed for random number generation to ensure exact reproducibility. |
| `re` | `"R1"` (`"R1"`, `"R2"`) | Random-effects structure as defined in the manuscript. |
| `sigma2_b` | `1` | Variance of the random effects ($\sigma_b^2$). |
| `sigma` | `1` | Error standard deviation ($\sigma$). |
| `rho_b` | `0.3` | Correlation parameter among random effects ($\rho_b$). |
| `signal` | `"S2"` (`"S1"`, `"S2"`) | True active predictor sparsity pattern:<br>• `"S1"`: active predictors at positions $1:5$<br>• `"S2"`: active predictors at positions $\{5, 10, 50, 100, 200\}$ |
| `kappa` | `0.2` (`0`, `0.2`, etc.) | Signal decay parameter ($B = n^{-\kappa}$):<br>• `0`: signal magnitude remains constant with sample size<br>• `0.2`: signal magnitude decreases with sample size |
| `sigx` | `"CS"` (`"I"`, `"CS"`, `"AR1"`) | Covariance structure of fixed-effect predictors $X$: Independent (`"I"`), Compound Symmetry (`"CS"`), or First-order Autoregressive (`"AR1"`). |
| `rho_x` | `0.5` | Off-diagonal predictor correlation parameter ($\rho_x$) for CS/AR1 structures. |
| `m` | `10` | Total number of subjects/clusters ($m$). |
| `ni` | `rep(c(10, 15, 20, 25, 30), 2)` | Vector or scalar specifying sample sizes per cluster ($n_i$). |
| `dist` | `"gauss"` (`"gauss"`, `"heavy"`, `"skew"`, `"RE_t5"`) | Distributional specification for errors and random effects:<br>• `"gauss"`: Normal errors + Normal random effects<br>• `"heavy"`: Heavy-tailed $t_5$ errors + Normal random effects<br>• `"skew"`: Skew-normal errors + Normal random effects<br>• `"RE_t5"`: Normal errors + Multivariate-$t_5$ random effects |
| `sizes` | `seq(10, 200, 20)` | Grid of retained model sizes ($d$) evaluated in variable screening. |


Contamination Schemes (C0)--(C5) are as specified in the manuscript (they can be changed under the settings `CONTAM`): 

* **$C0$**: Clean baseline data (0% contamination).
* **$C1$**: Casewise additive contamination on the response vector $y$ ($5\%$ and $10\%$).
* **$C2$**: Casewise additive contamination on the first active fixed-effect predictors ($5\%$ and $10\%$).
* **$C3$**: High-leverage outlier contamination on the first inactive fixed-effect predictors ($5\%$ and $10\%$).
* **$C4$**: High-leverage outlier contamination on random-effect design matrix $Z$ ($5\%$ and $10\%$).
* **$C5$**: Cluster-level block contamination affecting both response $y$ and predictors $X$ ($20\%$ clusters).

---



### Instructions for Execution

First, install the additionally required pakages for its own parallel execution, data generation,
reproducible parallel RNG streams, skew-normal error generation, and single-threaded BLAS timing, respectively:
```r
install.packages(c("foreach", "doParallel"))
install.packages(c("doRNG", "sn", "RhpcBLASctl"))
```

Also, install the package for data generation:
```r
install.packages("devtools")
devtools::install_github("abhik-stat/fsimR")
```

<br>

#### Option A: Running interactively in `R` / `RStudio`

1. Open the R session from the package root directory.
2. Edit configuration options in `simulation_script.R` if desired.
3. Execute the script:
```R
source("simulation_script.R")

```



#### Option B: Running from terminal / command line

```bash
# Execute in background with log output
Rscript simulation_script.R

```

---

### Generated Results and Output Structure

When `simulation_script.R` runs, it dynamically creates output directories and populates them as follows:

```text
.
├── <TAG>_<TIMESTAMP>.log         # High-precision execution log with timestamps
├── <TAG>_summary.csv             # Consolidated summary metrics (TPR, EmpSSP, FDP, MinMS)
├── Results (Raw)/                # Per-scheme raw replication data
│   ├── <TAG>_C0_prop00.csv
│   ├── <TAG>_C1_prop05.csv
│   ├── <TAG>_proxy_errors.csv   # Operator-norm error summaries for proxy estimators
│   └── ...
└── Figures/                      # Generated visualization plots
    ├── <TAG>_MinMS_boxplots.pdf  # Minimum model size boxplots per contamination cell
    ├── <TAG>_runtime_boxplots.pdf# Computational runtime comparisons
    └── <TAG>_TPRcurve.pdf        # Publication-ready TPR curves across retained model sizes (d)

```

*(Note: `<TAG>` is a auto-generated string summarizing the configuration, e.g., `S2_R1_CS(0.5)_RE(1_0.3)_m10_ni10-15-20-25-30_RE_t5_p1000_Bk2`)*


---



<br>
<br>

## Citation

If you use this code, please cite:

```bibtex
@article{ghosh_thoresen_dpdsisp,
  title   = {Robust and Scalable Sure Screening of Fixed Effects in Ultrahigh-dimensional
             Linear Mixed Models},
  author  = {Ghosh, Abhik and Thoresen, Magne},
  journal = {Statistics and Computing},
  year    = {2026},
  note    = {Under revision; update with volume/pages/DOI upon acceptance}
}
```

## Contact

For questions regarding the code or the associated paper, 
please contact Dr. Abhik Ghosh at abhik.ghosh.stat@gmail.com.

Bug reports, suggestions, and pull requests are welcome. 

