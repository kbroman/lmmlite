context("Standard errors on heritability are OK")

test_that("additionally computing SE for h2 does nothing surprising", {

    # expected results
    load("expected_lmm_results.RData") # expected_lmm_r and expected_lmm_m

    tol <- 1e-8 # tolerance for convergence

    # load test data
    data(recla)
    k <- recla$kinship
    y <- recla$pheno[,1]
    X <- recla$covar
    e <- eigen_rotation(k, y, X, use_cpp=TRUE)
  
    # Fit model for a given phenotype with and without standard errors
    m1 <- fitLMM(e$Kva, e$y, e$X, tol=tol, use_cpp=FALSE)
    m2 <- fitLMM(e$Kva, e$y, e$X, tol=tol, use_cpp=FALSE, compute_se = TRUE)
    
    # Verify legitmmate SE
    expect_gt(attr(m2$hsq, "se"), 0)
    
    # Verify equivalent output one removing the SE
    m2$hsq <- as.numeric(m2$hsq)
    expect_equal(m1, m2)
})
