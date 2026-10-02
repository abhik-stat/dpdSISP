# dpdSISP
Robust Sure Screening for Ultrahigh-Dimensional Linear Mixed Models using DPD-SISP

**DPD-SISP** is a minimum density power divergence (DPD) marginal screening procedure for
fixed effects in linear mixed models. It first *whitens* the data using a proxy for the
random-effects covariance structure, then ranks candidate covariates by a robust marginal
association statistic -- designed to remain reliable when a fraction of observations are
contaminated, and to scale to tens of thousands of candidates.

Code accompanying:

> Ghosh, A. and Thoresen, M. **Robust and Scalable Sure Screening of Fixed Effects in
> Ultrahigh-dimensional Linear Mixed Models.** *Statistics and Computing* (under revision).
