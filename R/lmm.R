# simple port of pyLMM to R

#' eigen decomposition + rotation
#'
#' Do eigen decomposition of kinship matrix and rotate \code{X} and
#' \code{y} by that, i.e., pre-multiply by the transpose of the matrix
#' of eigenvectors. If \code{Kva} and \code{Kve_t} provided, just do
#' the "rotation".
#'
#' @param K Kinship matrix
#' @param y Phenotypes
#' @param X Numeric matrix with covariates. If missing, use a column
#' of 1's (for intercept).
#' @param Kva Eigenvalues of \code{K} (optional)
#' @param Kve_t = transposed eigenvectors of K (optional)
#'
#' @export
#' @return List containing \code{Kva}, \code{Kve_t} and rotated
#' \code{y} and \code{X}.
#'
#' @examples
#' data(recla)
#' e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar)
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


#' Get MLEs for coefficients and variance
#'
#' For a fixed value for \code{hsq}, the heritability, calculate the
#' corresponding maximum likelihood estimates of \code{beta} and
#' \code{sigsq}, with the latter being the total variance,
#' \code{sigsq_g + sigsq_e}.
#'
#' @param hsq heritability
#' @param Kva eigenvalues of K (calculated by \code{\link{eigen_rotation}})
#' @param y rotated phenotypes (calculated by \code{\link{eigen_rotation}})
#' @param X rotated covariate matrix (calculated by \code{\link{eigen_rotation}})
#' @param reml If TRUE, use REML; otherwise use ordinary maximum likelihood.
#'
#' @export
#' @return list containing \code{beta} and \code{sigsq}, with residual
#' sum of squares and (if \code{reml=TRUE}, \code{log det (XSX)}) as
#' attributes.
#'
#' @examples
#' data(recla)
#' e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar)
#' ml <- getMLsoln(0.5, e$Kva, e$y, e$X)
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

    # calculate a bunch of matrices
    XSt = t(X * S)     # (XS)'
    ySt = t(y * S)     # (yS)'
    XSX = XSt %*% X    # (XS)'X
    XSy = XSt %*% y    # (XS)'y
    ySy = ySt %*% y    # (yS)'y

    # estimate of beta, by weighted LS
    e <- eigen(XSX)
    evals <- e$values
    evecs <- t(e$vectors)
    if(reml) logdetXSX <- sum(log(evals))
    beta <- t(evecs/evals) %*% evecs %*% XSy

    rss = ySy - t(XSy) %*% beta

    # estimate of sigma^2 (total variance = sigma_g^2 + sigma_e^2)
    sigsq <- rss / ifelse(reml, n - p, n)

    # return value
    result <- list(beta=beta, sigsq=sigsq)
    attr(result, "rss") <- rss
    if(reml) attr(result, "logdetXSX") <- logdetXSX

    result
}

#' Calculate log likelihood for a given heritability
#'
#' Calculate the log likelihood for a given value of the heritability, \code{hsq}.
#'
#' @param hsq heritability
#' @param Kva eigenvalues of K (calculated by \code{\link{eigen_rotation}})
#' @param y rotated phenotypes (calculated by \code{\link{eigen_rotation}})
#' @param X rotated covariate matrix (calculated by \code{\link{eigen_rotation}})
#' @param reml If TRUE, use REML; otherwise use ordinary maximum likelihood.
#'
#' @export
#' @return The log likelihood value, with the corresponding estimates
#' of \code{beta} and \code{sigsq} included as attributes.
#'
#' @examples
#' data(recla)
#' e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar)
#' loglik <- calcLL(0.5, e$Kva, e$y, e$X)
#' many_loglik <- calcLL(seq(0, 1, by=0.1), e$Kva, e$y, e$X)
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

    # calculate log likelihood
    rss <- attr(MLsoln, "rss")
    LL <- -0.5*(sum(log(hsq*Kva + 1-hsq)) + n*log(rss))

    if(reml) { # note that default is determinant() gives log det
        logdetXpX <- attr(X, "logdetXpX")
        if(is.null(logdetXpX)) { # need to calculate it
            XpX <- t(X) %*% X
            logdetXpX <- sum(log(eigen(XpX)$values))
        }

        logdetXSX <- attr(MLsoln, "logdetXSX")
        LL <- LL + 0.5 * (p*log(sigsq) + logdetXpX - logdetXSX)
    }

    attr(LL, "beta") <- beta
    attr(LL, "sigsq") <- sigsq
    LL
}

#' Fit a linear mixed model
#'
#' Fit a linear mixed model of the form y = Xb + e where e follows a
#' multivariate normal distribution with mean 0 and variance matrix
#' \code{sigsq_g K + sigsq_e I}, where \code{K} is a known kniship
#' matrix and \code{I} is the identity matrix.
#'
#' @param Kva Eigenvalues of K (calculated by \code{\link{eigen_rotation}})
#' @param y Rotated phenotypes (calculated by \code{\link{eigen_rotation}})
#' @param X Rotated covariate matrix (calculated by \code{\link{eigen_rotation}})
#' @param reml If TRUE, use REML; otherwise use ordinary maximum likelihood.
#' @param tol Tolerance for convergence
#'
#' @export
#' @return List containing estimates of \code{beta}, \code{sigsq},
#' \code{hsq}, \code{sigsq_g}, and \code{sigsq_e}, as well as the log
#' likelihood (\code{loglik}).
#'
#' @examples
#' data(recla)
#' e <- eigen_rotation(recla$kinship, recla$pheno[,1], recla$covar)
#' result <- fitLMM(e$Kva, e$y, e$X)
fitLMM <-
    function(Kva, y, X, reml=TRUE, tol=1e-4)
{
    n <- length(Kva)
    if(!is.matrix(X)) X <- as.matrix(X)
    stopifnot(nrow(X) == n)
    if(!is.matrix(y)) y <- as.matrix(y)
    stopifnot(nrow(y) == n)

    # calculate log determinant of X'X matrix, so it's only done once
    if(reml) {
        XpX <- t(X) %*% X
        attr(X, "logdetXpX") <- sum(log( eigen(XpX)$values ))
    }

    # maximize log likelihood
    out <- stats::optimize(calcLL, c(0, 1), Kva=Kva, y=y, X=X, reml=reml,
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
