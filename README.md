# dpdSISP: Robust Sure Screening for Linear Mixed Models

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

---


## General Function Library

Every screening procedures used in our study (Ghosh and Thoresen, 2026) is defined 
as a generic functions, along with necessary background and heper function, in `R/`. 


### Requirements

`R` (>= 4.1 recommended).

```r
install.packages(c("MASS", "mvtnorm", "lme4", "glmnet", "Matrix"))
```


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


<br>

---

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

* (C0): Clean baseline data (0% contamination).
* (C1): Casewise additive contamination on the response vector $y$ ($5\%$ and $10\%$).
* (C2): Casewise additive contamination on the first active fixed-effect predictors ($5\%$ and $10\%$).
* (C3): High-leverage outlier contamination on the first inactive fixed-effect predictors ($5\%$ and $10\%$).
* (C4): High-leverage outlier contamination on random-effect design matrix $Z$ ($5\%$ and $10\%$).
* (C5): Cluster-level block contamination affecting both response $y$ and predictors $X$ ($20\%$ clusters).

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


<br>

---

---


## Real Data Application: Reproducing ADNI-2 analyses from the paper 

This section provides instructions for reproducing the real data application results 
presented in the manuscript Ghosh and Thoresen (2026) using the master script 
`Data_analysis_script.R`. This script executes the complete analytical workflow 
integrating longitudinal cognitive scores (MMSE) and high-dimensional blood gene expression profiles 
from the Alzheimer's Disease Neuroimaging Initiative Phase 2(ADNI-2) cohort, 
performing required preprocessing and comparing various robust screening procedures 
and stability selection under high-dimensional LMM.

<br>

### Data Availability & Access Statement

**Important Note on Data Privacy:** 
The datasets utilized in this study contain restricted human subject health information 
from the Alzheimer's Disease Neuroimaging Initiative (ADNI) 
and cannot be shared or redistributed directly with this repository.

To run the provided R scripts, users must independently obtain data access permissions via teh following steps:

1. Register at [adni.loni.usc.edu](https://adni.loni.usc.edu). Request the `ADNIMERGE2` R
   package (Study Files -> Study Info -> Data & Database) and install it locally:
   `install.packages("path/to/ADNIMERGE2_<version>.tar.gz", repos = NULL, type = "source")`
   ([documentation](https://atri-biostats.github.io/ADNIMERGE2/)).
2. Separately download the gene-expression profile (Genetic Data -> "Gene Expression" study on
   [ida.loni.usc.edu](https://ida.loni.usc.edu)) and its accompanying probe manifest.
   Downlaod these two files, having names `ADNI_gene_expression_profile.csv` and `gene_probe_manifest.tsv`, respectively, 
   from the ADNI genomic data section and place them in the `rawData/` folder.

Once these files and packages are acquired, the supplied pipeline code (`Data_analysis_script.R`) can be executed directly to reproduce the analysis.


---

### Dataset Overview & Requirements

The analysis integrates two primary ADNI data sources:

1. **Longitudinal Clinical & Cognitive Data (`ADNIMERGE2` R Package):**
   - **Response ($y$):** Mini-Mental State Examination (MMSE) scores (`MMSCORE`) measured repeatedly across visits.
   - **Conditioning Covariates ($X_C$):** Baseline age, gender, education level (`EDUC`), baseline diagnosis (`DX`), and APOE $\varepsilon4$ allele count (`apoe4`).
   - **Time Trajectory ($Z$):** Visit time measured in years from baseline ($Z_0 = 1$ for random intercept, $Z_1 = \text{Time}$ for random slope).

2. **Gene Expression Profiles (`ADNI_gene_expression_profile.csv`):**
   - Whole-blood transcriptomic profiling using Affymetrix Human Genome U219 arrays (49,386 probe sets).
   - Matched against sample manifest metadata (`gene_probe_manifest.tsv`).

---

### Data Assembly & Preprocessing Pipeline

The assembly module executes the following sequential steps:

1. **Inclusion Criteria:**
   - Restricts cohort to the target study (`ADNI2`).
   - Filters participants with at least $3$ longitudinal visits ($\ge 3$ MMSE records).
   - Retains complete cases with respect to required conditioning covariates.

2. **Transcriptomic Quality Control & Batch Correction:**
   - **Batch Correction:** Applies Empirical Bayes batch correction (`sva::ComBat`) across Affymetrix plate identifiers.
   - **Control Filtering:** Removes Affymetrix internal control probes (`AFFX-*`).
   - **Annotation Matching:** Drops unmapped probes lacking official HGNC gene symbols.
   - **Variance Filtering:** Filters out non-informative probes falling below a specified empirical variance quantile.

3. **Sample Alignment & Matrix Construction:**
   - Aligns baseline gene expression arrays ($p$ probes) with repeated visit records ($n$ total visits across $m$ unique participants).
   - Generates both probe-level ($X_p$) and gene-level aggregated matrices ($X_{\text{gene}}$, retaining the probe with maximum interquartile range per gene).

---

### High-Dimensional Screening & Stability Selection

**Screening Methods**: The pipeline evaluates screening models by conditioning on clinical covariates ($X_C$) while screening transcriptomic predictors:

| Method Code | Type | Description |
| :--- | :--- | :--- |
| `MLE` / `REML` | Benchmark | Standard ML / REML estimator based screening under marginal LMMs|
| `TPCc` | Benchmark | Two-stage Partial Correlation Conditioning |
| `cv-P` ($\alpha$) | DPD (Proposed) | DPD-SISP using cross-validated proxy estimators ($\alpha \in \{0.1, 0.3, 0.5\}$) |
| `I0-P` ($\alpha$) | DPD (Proposed) | DPD-SISP using null intercept-only proxy estimators ($\alpha \in \{0.1, 0.3, 0.5\}$) |
| `no-P` ($\alpha$) | DPD (Baseline) | DPD-SIS ignoring random effects structure, i.e., DPD-SISP with null proxy ($\alpha \in \{0.1, 0.3, 0.5\}$) |

**Screen Size Conventions**: 
Screening is performed at target size $d = \lfloor n / \log n \rfloor$ (record-count convention).
Other applicable option is $d = \lfloor m / \log m \rfloor$ (subject-count convention).

**Complementary Pair Stability Selection (CPSS)**: 
To control false discoveries under high-dimensional noise and potential outliers:
- Evaluates selection frequencies across $B = 25$ subsampling splits ($2B = 50$ half-runs).
- Aggregates empirical selection probabilities $\hat{\pi}_j$ across methods to identify highly stable biomarker signatures.

---

### Instructions for Execution

First, install the additionally required packages for batch correction, gene expression preprocessing, parallel execution, and reproducible RNG streams:
```r
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager")

BiocManager::install("sva")                          # only if batch correction is turned on
install.packages(c("VennDiagram", "UpSetR"))         # used by summary_results.R's diagrams

install.packages(c("foreach", "doParallel", "doRNG"))
```

(Note: Ensure `ADNIMERGE2` is installed from the official ADNI portal as described above).

<br>



#### Option A: Running interactively in `R` / `RStudio`

1. Open the R session from the package root directory.
2. Edit configuration options in `Data_analysis_script.R`, if desired.
3. Execute the script:

```R
source("Data_analysis_script.R")
```

#### Option B: Running from terminal / command line

```bash
# Execute in background with log output
Rscript Data_analysis_script.R
```


### Generated Output Files

* `DerivedData_ADNI2_MMSCORE_...rds`: Cleaned long-format analysis dataset and expression matrices.
* `DerivedData_ADNI2_MMSCORE_..._screeningResults.rds`: Full-sample marginal screening ranks and top gene lists.
* `DerivedData_ADNI2_MMSCORE_..._cpss_results.rds`: Stability selection probabilities and selection frequency tables.
* `DerivedData_ADNI2_MMSCORE_..._pipeline_log.txt`: Timestamped execution log file.


<br>

---

---


### Citation

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

### License

All code and scripts in this repository are released under the [MIT
License](https://opensource.org/licenses/MIT),
which permits their reuse, modification, and distribution for academic or 
commercial purposes, provided that the original copyright notice and permission notice are included.

### Contact

For questions regarding the code or the associated paper, 
please contact **Dr. Abhik Ghosh** at *abhik.ghosh.stat@gmail.com*.

Bug reports, suggestions, and pull requests are welcome. 
