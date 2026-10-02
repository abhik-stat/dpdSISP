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
<br>



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

