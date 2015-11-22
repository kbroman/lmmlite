context("lmm in c++")

test_that("X'X works", {

    data(recla)

    expected <- t(recla$covar) %*% recla$covar
    dimnames(expected) <- NULL

    expect_equal(lmmlite:::R_calc_xpx(recla$covar), expected)

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

    # REML
    outR <- getMLsoln(0.5, e$Kva, e$y, e$X, TRUE)
    outcpp <- R_getMLsoln(0.5, e$Kva, e$y, e$X, TRUE)

    expect_equal(outcpp$beta, as.numeric(outR$beta))
    expect_equal(outcpp$sigsq, as.numeric(outR$sigsq))
    expect_equal(outcpp$logdetXSX, attr(outR, "logdetXSX"))
    expect_equal(outcpp$rss, as.numeric(attr(outR, "rss")))

    # real expected values
    expected <- list(sigsq=276158.504472848,
                     rss=71525052.6584677,
                     logdetXSX=6.96036161615051,
                     beta=c(1413.98690713607, -28.1793101783314))
    expect_equal(outcpp, expected)

})
