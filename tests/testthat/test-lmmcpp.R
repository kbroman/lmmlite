context("lmm in c++")

test_that("X'X works", {

    data(recla)

    expected <- t(recla$covar) %*% recla$covar
    dimnames(expected) <- NULL

    expect_equal(R_calc_xpx(recla$covar), expected)

})

test_that("eigen decomp works", {

    data(recla)

    e <- R_eigen_decomp(recla$kinship)

    eR <- eigen(recla$kinship)

    # eigenvalues the same as from R (but different order)
    expect_equal(sort(e$values), sort(eR$values))

    # matrix multiply and get original
    expected <- recla$kinship
    dimnames(expected) <- NULL
    expect_equal(t(e$vectors) %*% diag(e$values) %*% e$vectors, expected)

    # inverse
    expect_equal(t(e$vectors) %*% diag(1/e$values) %*% e$vectors, solve(expected))

})


test_that("getMLsoln works", {

    data(recla)
    e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar)

    # all that reml=TRUE does is calculate logdetXSX
    outR <- getMLsoln(0.5, e$Kva, e$y, e$X, TRUE)
    outcpp <- R_getMLsoln(0.5, e$Kva, e$y, e$X, TRUE)

    expect_equal(outcpp$beta, as.numeric(outR$beta))
    expect_equal(outcpp$sigmasq, as.numeric(outR$sigmasq))
    expect_equal(outcpp$logdetXSX, attr(outR, "logdetXSX"))
    expect_equal(outcpp$rss, as.numeric(attr(outR, "rss")))

    # real expected values
    expected <- list(sigmasq=276158.504472848,
                     rss=71525052.6584677,
                     logdetXSX=6.96036161615051,
                     beta=c(1413.98690713607, -28.1793101783314))
    expect_equal(outcpp, expected)

})


test_that("calcLL works", {

    data(recla)
    e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar)
    logdetXpX <- determinant(R_calc_xpx(e$X))$modulus

    #ML
    outR <- calcLL(0.5, e$Kva, e$y, e$X, FALSE)
    outcpp <- R_calcLL(0.5, e$Kva, e$y, e$X, FALSE)
    outcpp2 <- R_calcLL(0.5, e$Kva, e$y, e$X, FALSE, logdetXpX)

    expect_equal(outcpp$loglik, as.numeric(outR))
    expect_equal(outcpp$beta, as.numeric(attr(outR, "beta")))
    expect_equal(outcpp$sigmasq, as.numeric(attr(outR, "sigmasq")))

    expect_equal(outcpp2$loglik, as.numeric(outR))
    expect_equal(outcpp2$beta, as.numeric(attr(outR, "beta")))
    expect_equal(outcpp2$sigmasq, as.numeric(attr(outR, "sigmasq")))

    # real expected values
    expected <- list(loglik=-2349.20302794088,
                     sigmasq=276158.504472848,
                     beta=c(1413.98690713607, -28.1793101783314))
    expect_equal(outcpp, expected)
    expect_equal(outcpp2, expected)

    # REML
    outR <- calcLL(0.5, e$Kva, e$y, e$X, TRUE)
    outcpp <- R_calcLL(0.5, e$Kva, e$y, e$X, TRUE)
    outcpp2 <- R_calcLL(0.5, e$Kva, e$y, e$X, TRUE, logdetXpX)

    expect_equal(outcpp$loglik, as.numeric(outR))
    expect_equal(outcpp$beta, as.numeric(attr(outR, "beta")))
    expect_equal(outcpp$sigmasq, as.numeric(attr(outR, "sigmasq")))

    expect_equal(outcpp2$loglik, as.numeric(outR))
    expect_equal(outcpp2$beta, as.numeric(attr(outR, "beta")))
    expect_equal(outcpp2$sigmasq, as.numeric(attr(outR, "sigmasq")))

    # real expected values
    expected <- list(loglik=-2333.44582306843,
                     sigmasq=276158.504472848,
                     beta=c(1413.98690713607, -28.1793101783314))
    expect_equal(outcpp, expected)
    expect_equal(outcpp2, expected)

})
