## [R/lmmlite](https://kbroman.org/lmmlite) - R port of [pylmm](https://github.com/nickFurlotte/pylmm)

[![Build Status](https://travis-ci.org/kbroman/lmmlite.svg?branch=master)](https://travis-ci.org/kbroman/lmmlite)

Karl Broman (following the code in
[Nick Furlotte](http://whatmind.com)'s [pylmm](https://github.com/nickFurlotte/pylmm))

In order to better understand [pylmm](https://github.com/nickFurlotte/pylmm),
for linear mixed models for genome-wide association studies
(GWAS), I wrote an R version of the algorithm.

---

### Vignettes

- [User guide](https://kbroman.org/lmmlite/assets/lmmlite.html)
  [[source](https://github.com/kbroman/lmmlite/blob/master/vignettes/lmmlite.Rmd)]

- [Results compared](https://kbroman.org/lmmlite/assets/compare2pylmm.html)
  to [pylmm](https://github.com/nickFurlotte/pylmm) (and also the R package
  [regress](https://cran.r-project.org/web/packages/regress/)).
  [[source](https://github.com/kbroman/lmmlite/blob/gh-pages/assets/compare2pylmm.Rmd)]

- [Performance compared](https://kbroman.org/lmmlite/assets/performance.html)
  to pylmm.
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
[GNU Affero GPL](https://www.gnu.org/licenses/why-affero-gpl.html).

The code was developed following study of [Nick Furlotte](http://whatmind.com)'s
[pylmm](https://github.com/nickFurlotte/pylmm) code.
