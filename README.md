# dpdSISP
Robust Sure Screening for Ultrahigh-Dimensional Linear Mixed Models using DPD-SISP

**DPD-SISP** is a minimum density power divergence (DPD) estimator based marginal screening procedure for
fixed-effects in linear mixed models. It first *whitens* the data using a proxy for the unknown
random-effects covariance structure, then ranks candidate covariates by a robust marginal
association statistic, designed to remain reliable when a fraction of observations are
contaminated, and to scale to tens of thousands of candidates.

Code accompanying:

> Ghosh, A. and Thoresen, M. **Robust and Scalable Sure Screening of Fixed Effects in
> Ultrahigh-dimensional Linear Mixed Models,** which is presently under review in *Statistics and Computing*.



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

### Contact

For questions regarding the code or the associated paper, 
please contact Dr. Abhik Ghosh at abhik.ghosh.stat@gmail.com.

Bug reports, suggestions, and pull requests are welcome. 

