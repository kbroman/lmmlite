## [R/lmmlite](https://kbroman.org/lmmlite) - R port of [pylmm](https://github.com/nickFurlotte/pylmm)

[![R-CMD-check](https://github.com/kbroman/lmmlite/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/kbroman/lmmlite/actions/workflows/R-CMD-check.yaml)
[![zenodo DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.5149506.svg)](https://doi.org/10.5281/zenodo.5149506)

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
  [regress](https://cran.r-project.org/package=regress)).
  [[source](https://github.com/kbroman/lmmlite/blob/gh-pages/assets/compare2pylmm.Rmd)]

- [Performance compared](https://kbroman.org/lmmlite/assets/performance.html)
  to pylmm.
  [[source](https://github.com/kbroman/lmmlite/blob/gh-pages/assets/performance.Rmd)]


---

### Installation

You can install R/lmmlite from
[GitHub](https://github.com/kbroman/lmmlite).

You first need to install the
[remotes](https://remotes.r-lib.org) pakcage.

    install.packages("remotes")

Then use `install_github()` to install R/lmmlite.

    library(remotes)
    install_github("kbroman/lmmlite")

The [Rcpp](https://github.com/RcppCore/Rcpp) and
[RcppEigen](https://github.com/RcppCore/RcppEigen) packages
will also be installed.

---

### License

[R/lmmlite](https://github.com/kbroman/lmmlite) is released under the
[GNU Affero GPL](https://www.gnu.org/licenses/why-affero-gpl.html).

The code was developed following study of [Nick Furlotte](http://whatmind.com)'s
[pylmm](https://github.com/nickFurlotte/pylmm) code.
