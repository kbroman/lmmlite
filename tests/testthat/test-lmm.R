context("lmm regression test")

test_that("fitLMM works", {

    # expected results
    load("expected_lmm_results.RData") # expected_lmm_r and expected_lmm_m

    tol <- 1e-8 # tolerance for convergence

    # load test data
    data(recla)
    k <- recla$kinship
    y <- recla$pheno
    X <- recla$covar

    # scale phenotypes to have sd 1
    y <- t( t(y) / apply(y, 2, sd, na.rm=TRUE) )

    # analyses all phenotypes by reml
    lmm_all_r <- lapply(1:ncol(y), function(i) {
        thisy <- y[,i,drop=FALSE]
        omit <- is.na(thisy)
        thisy <- thisy[!omit,,drop=FALSE]
        thisX <- X[!omit,,drop=FALSE]
        thisk <- k[!omit,!omit]
        e <- eigen_rotation(thisk, thisy, thisX, use_cpp=FALSE)
        fitLMM(e$Kva, e$y, e$X, tol=tol, use_cpp=FALSE)})

    # combine results
    tab_lmm_r <- t(vapply(lmm_all_r, function(a) c(sigmasq_g=a$sigmasq_g,
                                                   sigmasq_e=a$sigmasq_e,
                                                   hsq=a$hsq,
                                                   beta_int=a$beta[1],
                                                   beta_sex=a$beta[2],
                                                   loglik=a$loglik),
                          rep(0, 6)))

    # compare to expectedd
    expect_equal(tab_lmm_r, expected_lmm_r)


    # analyses all phenotypes by ML
    lmm_all_m <- lapply(1:ncol(y), function(i) {
        thisy <- y[,i,drop=FALSE]
        omit <- is.na(thisy)
        thisy <- thisy[!omit,]
        thisX <- X[!omit,]
        thisk <- k[!omit,!omit]
        e <- eigen_rotation(thisk, thisy, thisX, use_cpp=FALSE)
        fitLMM(e$Kva, e$y, e$X, tol=tol, reml=FALSE, use_cpp=FALSE)})

    # combine results
    tab_lmm_m <- t(vapply(lmm_all_m, function(a) c(sigmasq_g=a$sigmasq_g,
                                                   sigmasq_e=a$sigmasq_e,
                                                   hsq=a$hsq,
                                                   beta_int=a$beta[1],
                                                   beta_sex=a$beta[2],
                                                   loglik=a$loglik),
                          rep(0, 6)))

    # compare to expected
    expect_equal(tab_lmm_m, expected_lmm_m)
})
