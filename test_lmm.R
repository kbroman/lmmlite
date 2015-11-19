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
expect_equal(out1m$sigma, grab_sigs(lmm1m), tolerance=0.00001)
