# code to test lmm.R
tol <- 1e-8 # tolerance for convergence

source("lmm.R")

source("load_test_data.R")
# scale phenotypes to have sd 1
y <- t( t(y) / apply(y, 2, sd, na.rm=TRUE) )

# LMM by REML and ML with regress package
library(regress)
out1r <- regress(y[,1] ~ -1 + X, ~k, tol=tol)

# LMM with lmm.R
e <- eigen_rotation(k, y[,1], X)
lmm1r <- fitLMM(e$Kva, e$y, e$X, tol=tol)

# compare results
library(testthat)
expect_equal(out1r$sigma, c(k=lmm1r$sigsq_g, In=lmm1r$sigsq_e),
             tolerance=0.00001)
rownames(out1r$beta) <- gsub("^X", "", rownames(out1r$beta))
expect_equal(out1r$beta, lmm1r$beta, tolerance=0.0000001)

# do all phenotypes by reml
library(parallel)
lmm_all_r <- mclapply(1:ncol(y), function(i) {
    thisy <- y[,i,drop=FALSE]
    omit <- is.na(thisy)
    thisy <- thisy[!omit,,drop=FALSE]
    thisX <- X[!omit,,drop=FALSE]
    thisk <- k[!omit,!omit]
    e <- eigen_rotation(thisk, thisy, thisX)
    fitLMM(e$Kva, e$y, e$X, tol=tol)},
                     mc.cores=detectCores())

tab_lmm_r <- t(vapply(lmm_all_r, function(a) c(sigsq_g=a$sigsq_g,
                                               sigsq_e=a$sigsq_e,
                                               hsq=a$hsq,
                                               beta_int=a$beta[1],
                                               beta_sex=a$beta[2],
                                               loglik=a$loglik),
                      rep(0, 6)))



lmm_all_m <- mclapply(1:ncol(y), function(i) {
    thisy <- y[,i,drop=FALSE]
    omit <- is.na(thisy)
    thisy <- thisy[!omit,]
    thisX <- X[!omit,]
    thisk <- k[!omit,!omit]
    e <- eigen_rotation(thisk, thisy, thisX)
    fitLMM(e$Kva, e$y, e$X, tol=tol, reml=FALSE)},
                     mc.cores=detectCores())

tab_lmm_m <- t(vapply(lmm_all_m, function(a) c(sigsq_g=a$sigsq_g,
                                               sigsq_e=a$sigsq_e,
                                               hsq=a$hsq,
                                               beta_int=a$beta[1],
                                               beta_sex=a$beta[2],
                                               loglik=a$loglik),
                      rep(0, 6)))

source("expected_lmm.R")
expect_equal(tab_lmm_r, expected_lmm_r)
expect_equal(tab_lmm_m, expected_lmm_m)
