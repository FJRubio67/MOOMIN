# Moment-Objective Minimum-Discrepancy (MOOMIN) Prior

<div align="center">
<img width="400" alt="MOOMIN" src="https://github.com/user-attachments/assets/b2b65d4b-1723-4340-ae8a-5486aacc64e2" />
</div>

## Overview

This repository contains the R code used to reproduce all results in:

> Rubio, F.J. (2026). An objective non-local prior for skew-symmetric models.
> *Statistics & Probability Letters*, in press.
> [DOI: 10.1016/j.spl.2026.110730](https://doi.org/10.1016/j.spl.2026.110730) |
> [Preprint (arXiv)](https://arxiv.org/abs/2603.08285)

Skew-symmetric models extend symmetric distributions (such as the normal) by
introducing a shape parameter that controls skewness. However, standard objective
priors — such as Jeffreys' prior — tend to concentrate mass near the symmetric
special case, making it difficult to detect or estimate skewness from data.
The **MOOMIN prior** is an objective non-local prior for the shape parameter of
skew-symmetric models, constructed to assign negligible prior mass near symmetry
while remaining automatic (requiring no user-specified tuning). This facilitates
Bayesian model selection between symmetric and skew-symmetric models and improves
estimation of the shape parameter.

## Repository structure

```
MOOMIN/
├── routines/
│   └── routines_sn.R     # Core routines for real data application and simulations
├── Plots/                 # R code to reproduce all figures in the manuscript and appendix
└── Application/          # R script, R Markdown source, and rendered HTML for the
                          #   real data application
```

## Requirements

The following R packages are required to run the code:

```r
install.packages(c("sn", "MCMCpack", "coda", "ggplot2"))
```

> **Note:** exact package requirements can be verified by inspecting the
> `library()` calls at the top of `routines/routines_sn.R`.

## Tutorials and companion notes

- [Objective non-local priors for skew-symmetric models](https://rpubs.com/FJRubio/MOOMIN) — illustrative tutorial accompanying the paper
- [The effect of the shape parameter in skew-symmetric models](https://rpubs.com/FJRubio/DTVMINSS) — motivation via total variation distance
- [The effect of the shape parameter, Part II](https://rpubs.com/FJRubio/DiscMinSS) — discrepancy minimisation perspective
- [The effect of the shape parameter, Part III](https://rpubs.com/FJRubio/DivMinLC) — divergence minimisation and log-concavity

## Citation

If you use this code, please cite:

```bibtex
@article{rubio:2026,
  author  = {Rubio, F.J.},
  title   = {An objective non-local prior for skew-symmetric models},
  journal = {Statistics \& Probability Letters},
  year    = {2026},
  volume  = {235},
  pages   = {110730},
  doi     = {10.1016/j.spl.2026.110730}
}
```

## License

This repository is licensed under the [MIT License](LICENSE).
