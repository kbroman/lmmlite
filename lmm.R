# simple port of pyLMM to R

# do eigen decomposition of kinship matrix
# and rotate X and y by that [ie., pre-multiply by t(eigenvec)]
#
# If Kva and Kve_t provided, just do the "rotation"
#
# K = kinship matrix
# y = phenotypes
# X = covariate matrix
# Kva = eigenvalues of K (optional)
# Kve_t = transposed eigenvectors of K (optional)
eigen_rotation <-
    function(K, y, X, Kva, Kve_t)
{
    if(missing(Kva) || is.null(Kva) ||
       missing(Kve_t) || is.null(Kve_t)) {

        n <- nrow(K) # no. individuals
        stopifnot(ncol(K) == n) # square?

        # calculate eigen vals and vecs
        e <- eigen(K)
        Kva <- e$values
        Kve_t <- t(e$vectors)
    }
    else {
        n <- nrow(Kve_t)
        stopifnot(ncol(Kve_t) == n) # square?
        stopifnot(length(Kva) == n)
    }

    # more checks
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == n)

    if(missing(X) || is.null(X))
        X <- matrix(1, nrow=n)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == n)

    # rotation
    y <- Kve_t %*% y
    X <- Kve_t %*% X

    list(Kva=Kva, Kve_t=Kve_t, y=y, X=X)
}


# get MLEs for beta and sigsq
#
# hsq = heritability
# Kva = eigenvalues of K
# y = phenotypes (rotated)
# X = covariate matrix (rotated)
getMLsoln <-
    function(hsq, Kva, y, X, reml=TRUE)
{
    n <- length(Kva)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == n)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == n)
    p <- ncol(X)

    # diagonal matrix of weights
    S = 1/(hsq*Kva + 1-hsq)

    # X'S
    Xt = t(X * S)

    # X'SX
    XX = Xt %*% X

    # estimate of beta, by weighted LS
    # (use eigen decomposition of XX here, for later determinant?)
    beta = solve(XX, Xt %*% y)

    # resid
    resid = y - X %*% beta

    # RSS
    Q = sum(resid * S * resid)

    # estimate of sigma^2 (total variance = sigma_g^2 + sigma_e^2)
    sigsq <- Q / ifelse(reml, n - p, n)

    # return value
    result <- list(beta=beta, sigsq=sigsq)
    attr(result, "Q") <- Q
    attr(result, "XX") <- XX

    result
}

# calculate log likelihood of hsq
#
# hsq = heritability
# Kva = diagonal matrix of eigenvalues of K
# y = phenotypes (rotated)
# X = covariate matrix (rotated)
calcLL <-
    function(hsq, Kva, y, X, reml=TRUE)
{
    if(length(hsq) > 1)
        return(vapply(hsq, calcLL, 0, Kva, y, X, reml))

    n <- nrow(X)
    p <- ncol(X)

    # estimate beta and sigmasq
    MLsoln <- getMLsoln(hsq, Kva, y, X, reml=reml)
    beta <- MLsoln$beta
    sigsq <- MLsoln$sigsq
    Q <- attr(MLsoln, "Q")
    XX <- attr(MLsoln, "XX")

    # calculate log likelihood
    LL <- -0.5*(sum(log(hsq*Kva + 1-hsq)) + n*log(Q))

    if(reml) # note that default is determinant() gives log det
        LL <- LL + 0.5 * (p*log(sigsq) + determinant(t(X) %*% X)$modulus - determinant(XX)$modulus)

    attr(LL, "beta") <- beta
    attr(LL, "sigsq") <- sigsq
    LL
}

# fit LMM
#
# Kva = eigenvalues of K
# y = phenotypes (rotated)
# X = covariate matrix (rotated)
fitLMM <-
    function(Kva, y, X, reml=TRUE, tol=1e-4)
{
    n <- length(Kva)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == n)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == n)

    # maximize log likelihood
    out <- optimize(calcLL, c(0, 1), Kva=Kva, y=y, X=X, reml=reml,
                    maximum=TRUE, tol=tol)

    hsq <- out$maximum
    obj <- out$objective
    sigsq <- attr(obj, "sigsq")

    list(beta=attr(obj, "beta"),
         sigsq=sigsq, # total var
         hsq=hsq,
         sigsq_g= hsq*sigsq, # genetic variance
         sigsq_e = (1-hsq)*sigsq, # residual variance
         loglik = as.numeric(obj))
}
