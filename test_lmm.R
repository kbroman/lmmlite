# code to test lmm.R
tol <- 1e-8 # tolerance for convergence

source("lmm.R")

source("load_test_data.R")
# scale phenotypes to have sd 1
y <- t( t(y) / apply(y, 2, sd, na.rm=TRUE) )

# LMM by REML and ML with regress package
library(regress)
out1r <- regress(y[,1] ~ -1 + X, ~k, tol=tol)
out1m <- regress(y[,1] ~ -1 + X, ~k, kernel=0, tol=tol)

# LMM with lmm.R
e <- eigen_rotation(k, y[,1], X)
lmm1r <- fitLMM(e$Kva, e$y, e$X, tol=tol)
lmm1m <- fitLMM(e$Kva, e$y, e$X, reml=FALSE, tol=tol)

# grab sigma^2 values from lmm_result
grab_sigs <-
    function(lmm_result)
{
    c(k=lmm_result$sigsq_g,
      In=lmm_result$sigsq_e)
}


# compare results
library(testthat)
expect_equal(out1r$sigma, grab_sigs(lmm1r), tolerance=0.00001)
rownames(out1r$beta) <- gsub("^X", "", rownames(out1r$beta))
expect_equal(out1r$beta, lmm1r$beta, tolerance=0.0000001)

# do all phenotypes by reml
library(parallel)
lmmr_all <- mclapply(1:ncol(y), function(i) {
    thisy <- y[,i,drop=FALSE]
    omit <- is.na(thisy)
    thisy <- thisy[!omit,,drop=FALSE]
    thisX <- X[!omit,,drop=FALSE]
    thisk <- k[!omit,!omit]
    e <- eigen_rotation(thisk, thisy, thisX)
    fitLMM(e$Kva, e$y, e$X, tol=tol)},
                     mc.cores=parallel::detectCores())

regr_all <- mclapply(1:ncol(y), function(i) regress(y[,i] ~ X, ~k, tol=tol), mc.cores=parallel::detectCores())

sig_lmm <- vapply(lmmr_all, grab_sigs, c(0,0))
sig_regr <- vapply(regr_all, function(a) a$sigma, c(0,0))
plot(sig_lmm[1,], sig_regr[1,])
abline(0,1)

beta_lmm <- vapply(lmmr_all, function(a) a$beta, c(0,0))
beta_regr <- vapply(regr_all, function(a) a$beta, c(0,0))
plot((beta_lmm[1,] + beta_regr[1,])/2, beta_lmm[1,]- beta_regr[1,])
abline(0,1)




expect_equal(out1m$sigma, grab_sigs(lmm1m), tolerance=0.00001)
rownames(out1m$beta) <- gsub("^X", "", rownames(out1m$beta))
expect_equal(out1m$beta, lmm1m$beta, tolerance=0.0000001)

lmmm_all <- mclapply(1:ncol(y), function(i) {
    thisy <- y[,i,drop=FALSE]
    omit <- is.na(thisy)
    thisy <- thisy[!omit,]
    thisX <- X[!omit,]
    thisk <- k[!omit,!omit]
    e <- eigen_rotation(thisk, thisy, thisX)
    fitLMM(e$Kva, e$y, e$X, tol=tol, reml=FALSE)},
                     mc.cores=parallel::detectCores())
