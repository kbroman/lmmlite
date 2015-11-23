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
