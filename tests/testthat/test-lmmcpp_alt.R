context("lmm in c++ using plain R functions")

test_that("eigen_rotation works", {

    y <- recla$pheno[,1,drop=FALSE]
    X <- recla$covar
    dimnames(X) <- dimnames(y) <- NULL

    e <- eigen_rotation(recla$kinship, y, X, use_cpp=TRUE)
    expected <- eigen_rotation(recla$kinship, y, X, use_cpp=FALSE)

    # eigenvalues match
    expect_equal(sort(e$Kva), sort(expected$Kva))

    # back-rotation match
    expect_equal(t(e$Kve) %*% e$y, y)
    expect_equal(t(e$Kve) %*% e$X, X)

    # fit of LMM match
    out <- fitLMM(e$Kva, e$y, e$X, use_cpp=TRUE)
    outR <- fitLMM(expected$Kva, expected$y, expected$X, use_cpp=FALSE)

    expect_equal(out$loglik, as.numeric(outR$loglik))
    expect_equal(out$hsq, outR$hsq)
    expect_equal(out$beta, as.numeric(outR$beta))
    expect_equal(out$sigmasq, as.numeric(outR$sigmasq))

    # real expected values
    expected <- list(beta=c(1417.98330225743, -23.242180295058),
                     sigmasq=295340.52598977,
                     hsq=0.764284086070972,
                     sigmasq_g=225724.064011417,
                     sigmasq_e=69616.4619846284,
                     loglik=-2332.84011782658)
    expect_equal(out, expected)

})



test_that("getMLsoln works", {

    data(recla)
    e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar, use_cpp=TRUE)

    # all that reml=TRUE does is calculate logdetXSX
    outR <- getMLsoln(0.5, e$Kva, e$y, e$X, TRUE, use_cpp=FALSE)
    outcpp <- getMLsoln(0.5, e$Kva, e$y, e$X, TRUE, use_cpp=TRUE)

    expect_equal(outcpp$beta, as.numeric(outR$beta))
    expect_equal(outcpp$sigmasq, as.numeric(outR$sigmasq))
    expect_equal(attr(outcpp, "logdetXSX"), attr(outR, "logdetXSX"))
    expect_equal(attr(outcpp, "rss"), as.numeric(attr(outR, "rss")))

    # real expected values
    expected <- list(beta=c(1413.98690713607, -28.1793101783314),
                     sigmasq=276158.504472848)
    attr(expected, "logdetXSX") <- 6.96036161615051
    attr(expected, "rss") <- 71525052.6584677

    expect_equal(outcpp, expected)

})


test_that("calcLL works", {

    data(recla)
    e <- eigen_rotation(recla$kinship, recla$pheno[,1,drop=FALSE], recla$covar, use_cpp=TRUE)
    logdetXpX <- attr(e$X, "logdetXpX")

    #ML
    outR <- calcLL(0.5, e$Kva, e$y, e$X, FALSE, use_cpp=FALSE)
    outcpp <- calcLL(0.5, e$Kva, e$y, e$X, FALSE)

    expect_equal(as.numeric(outcpp), as.numeric(outR))
    expect_equal(attr(outcpp, "beta"), as.numeric(attr(outR, "beta")))
    expect_equal(attr(outcpp, "sigmasq"), as.numeric(attr(outR, "sigmasq")))

    # real expected values
    expected <- -2349.20302794088
    attr(expected, "beta") <- c(1413.98690713607, -28.1793101783314)
    attr(expected, "sigmasq") <- 276158.504472848
    expect_equal(outcpp, expected)

    # REML
    outR <- calcLL(0.5, e$Kva, e$y, e$X, TRUE, use_cpp=FALSE)
    outcpp <- calcLL(0.5, e$Kva, e$y, e$X, TRUE)

    expect_equal(as.numeric(outcpp), as.numeric(outR))
    expect_equal(attr(outcpp, "beta"), as.numeric(attr(outR, "beta")))
    expect_equal(attr(outcpp, "sigmasq"), as.numeric(attr(outR, "sigmasq")))

    # real expected values
    expected <- -2333.44582306843
    attr(expected, "beta") <- c(1413.98690713607, -28.1793101783314)
    attr(expected, "sigmasq") <- 276158.504472848

    expect_equal(outcpp, expected)

})

test_that("fitLMM works", {

    data(recla)
    e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar, use_cpp=TRUE)
    logdetXpX <- attr(e$X, "logdetXpX")

    #ML
    outR <- fitLMM(e$Kva, e$y, e$X, FALSE, use_cpp=FALSE)
    outcpp <- fitLMM(e$Kva, e$y, e$X, FALSE)

    expect_equal(outcpp$loglik, as.numeric(outR$loglik))
    expect_equal(outcpp$hsq, outR$hsq)
    expect_equal(outcpp$beta, as.numeric(outR$beta))
    expect_equal(outcpp$sigmasq, as.numeric(outR$sigmasq))

    # real expected values
    expected <- list(beta=c(1417.33690798302, -23.9850529710706),
                     sigmasq=291566.300541588,
                     hsq=0.720556143620898,
                     sigmasq_g=210089.889214325,
                     sigmasq_e=81476.4113477076,
                     loglik=-2348.79972912762)
    expect_equal(outcpp, expected)

    # REML
    outR <- fitLMM(e$Kva, e$y, e$X, TRUE, use_cpp=FALSE)
    outcpp <- fitLMM(e$Kva, e$y, e$X, TRUE)

    expect_equal(outcpp$loglik, as.numeric(outR$loglik))
    expect_equal(outcpp$hsq, outR$hsq)
    expect_equal(outcpp$beta, as.numeric(outR$beta))
    expect_equal(outcpp$sigmasq, as.numeric(outR$sigmasq))

    # real expected values
    expected <- list(beta=c(1417.98330225743, -23.242180295058),
                     sigmasq=295340.52598977,
                     hsq=0.764284086070972,
                     sigmasq_g=225724.064011417,
                     sigmasq_e=69616.4619846284,
                     loglik=-2332.84011782658)
    expect_equal(outcpp, expected)
})

test_that("fitLMM works for all phenotypes", {

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
        e <- eigen_rotation(thisk, thisy, thisX, use_cpp=TRUE)
        fitLMM(e$Kva, e$y, e$X, tol=tol, use_cpp=TRUE)})

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
        e <- eigen_rotation(thisk, thisy, thisX, use_cpp=TRUE)
        fitLMM(e$Kva, e$y, e$X, tol=tol, reml=FALSE, use_cpp=TRUE)})

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
