## [R/lmmlite](http://kbroman.org/lmmlite) - R port of [pylmm](https://github.com/nickFurlotte/pylmm)

[![Build Status](https://travis-ci.org/kbroman/lmmlite.svg?branch=master)](https://travis-ci.org/kbroman/lmmlite)

Karl Broman (following [Nick Furlotte](http://whatmind.com))

In order to better understand [pylmm](https://github.com/nickFurlotte/pylmm),
for linear mixed models for genome-wide association studies
(GWAS), I'm porting it to R.  The goal is learning, and ultimately
LMMs for QTL mapping.

Vignettes:

- [User guide](http://kbroman.org/lmmlite/assets/lmmlite.html)

- [Results compared](http://kbroman.org/lmmlite/assets/compare2pylmm.html)
  to [pylmm](https://github.com/nickFurlotte/pylmm) (and also the R package
  [regress](https://cran.r-project.org/web/packages/regress/)).

- [Performance compared](http://kbroman.org/lmmlite/assets/performance.html)
  to pylmm.


---

### Installation

You can install R/lmmlite from
[GitHub](https://github.com/kbroman/lmmlite).

You first need to install the
[devtools](https://github.com/hadley/devtools),
[Rcpp](https://github.com/RcppCore/Rcpp), and
[RcppEigen](https://github.com/RcppCore/RcppEigen) packages.

    install.packages(c("devtools", "Rcpp", "RcppEigen"))

Then use `devtools::install_github()` to install R/lmmlite.

    library(devtools)
    install_github("kbroman/lmmlite")

---

### License

[R/lmmlite](https://github.com/kbroman/lmmlite) is released under the
[GNU Affero GPL](https://www.gnu.org/licenses/why-affero-gpl.html),
because that's the license for
[pylmm](https://github.com/nickFurlotte/pylmm).
