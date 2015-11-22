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
