---
layout: page
title: R/lmmlite
description: R package for fitting linear mixed models for GWAS and QTL mapping
---

R/lmmlite is an R package for fitting linear mixed models for
genome-wide association studies (GWAS) and quantitative trait locus
(QTL) mapping.

I'm following the code in [pylmm](https://github.com/nickFurlotte/pylmm),
and really I'm developing this package mostly so that I can better understand pylmm.

The code is on [GitHub](https://github.com/kbroman/lmmlite).

---

### Vignettes

- [User guide](assets/lmmlite.html)
[[source](https://github.com/kbroman/lmmlite/blob/master/vignettes/lmmlite.Rmd)]

- [Comparison of results](assets/compare2pylmm.html) to those of
[pylmm](https://github.com/nickFurlotte/pylmm) and
[regress](https://cran.r-project.org/web/packages/regress/).
[[source](https://github.com/kbroman/lmmlite/blob/gh-pages/assets/compare2pylmm.Rmd)]

- [Performance comparison](assets/performance.html) with
[pylmm](https://github.com/nickFurlotte/pylmm).
[[source](https://github.com/kbroman/lmmlite/blob/gh-pages/assets/performance.Rmd)]

---

### Installation

You can install R/lmmlite from
[GitHub](https://github.com/kbroman/lmmlite).

You first need to install the
[devtools](https://github.com/hadley/devtools)
and [RcppEigen](https://github.com/RcppCore/RcppEigen) packages.
([Rcpp](https://github.com/RcppCore/Rcpp) will also be installed.)

    install.packages(c("devtools", "RcppEigen"))

Then use `devtools::install_github()` to install R/lmmlite.

    library(devtools)
    install_github("kbroman/lmmlite")

---

### License

[R/lmmlite](https://github.com/kbroman/lmmlite) is released under the
[GNU Affero GPL](https://www.gnu.org/licenses/why-affero-gpl.html),
because that's the license for
[pylmm](https://github.com/nickFurlotte/pylmm).