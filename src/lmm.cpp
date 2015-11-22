// linear mixed model via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include <math.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "lmm.h"

// calc X'X
MatrixXd calc_xpx(const MatrixXd& X)
{
    int n = X.cols();

    return MatrixXd(n,n).setZero().selfadjointView<Lower>()
        .rankUpdate(X.transpose());
}

// calc X'X (version to be called from R)
// [[Rcpp::export]]
NumericMatrix R_calc_xpx(const NumericMatrix& X)
{
    MatrixXd XX(as<Map<MatrixXd> >(X));
    return wrap(calc_xpx(XX));
}



// eigen decomposition
// returns eigenvalues and *transposed* eigenvectors
std::pair<VectorXd, MatrixXd> eigen_decomp(MatrixXd A)
{
    const SelfAdjointEigenSolver<MatrixXd> VLV(A);
    return std::make_pair(VLV.eigenvalues(), VLV.eigenvectors().transpose());
}

// eigen decomposition (version to be called from R)
// returns eigenvalues and *transposed* eigenvectors
// [[Rcpp::export]]
List R_eigen_decomp(NumericMatrix A)
{
    MatrixXd AA(as<Map<MatrixXd> >(A));
    std::pair<VectorXd,MatrixXd> result = eigen_decomp(AA);
    return List::create(Named("values") = result.first,
                        Named("vectors") = result.second);
}


// getMLsoln
// for fixed value of hsq, calculate MLEs of beta and sigmasq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = whether you'll be using REML (so need to calculate log det XSX)
VectorXd getMLsoln(double hsq, VectorXd Kva, VectorXd y,
                   MatrixXd X, bool reml=true)
{
    const int n = Kva.size();
    const int p = X.cols();

    // diagonal matrix of weights
    VectorXd S(n);
    for(int i=0; i<n; i++)
        S[i] = 1.0/(hsq*Kva[i] + 1.0-hsq);

    // calculate a bunch of matrices
    MatrixXd XSt = X.transpose() * S.asDiagonal();
    MatrixXd ySt(1,n);
    for(int i=0; i<n; i++) ySt(0,i) = y[i]*S[i];
    MatrixXd XSX = XSt * X;
    MatrixXd XSy = XSt * y;
    MatrixXd ySy = ySt * y;

    // estimate of beta, by weighted LS
    std::pair<VectorXd, MatrixXd>e = eigen_decomp(XSX);
    double logdetXSX=0.0;
    VectorXd inv_evals(p);
    for(int i=0; i<p; i++) {
        inv_evals[i] = 1.0/e.first[i];
        if(reml) logdetXSX += log(e.first[i]);
    }
    MatrixXd beta = e.second.transpose() * inv_evals.asDiagonal() * e.second * XSy;

    // residual sum of squares
    MatrixXd rss = ySy - XSy.transpose() * beta;

    // return value
    VectorXd result(p+3);
    result[0] = rss(0,0) / (double)(n-p); // sigma^2
    result[1] = rss(0,0);
    result[2] = logdetXSX; // determinant (if REML)
    for(int i=0; i<p; i++)
        result[i+3] = beta(i,0);

    return result;
}

// getMLsoln (version called from R)
// [[Rcpp::export]]
List R_getMLsoln(double hsq, NumericVector Kva, NumericVector y,
                 NumericMatrix X, bool reml=true)
{
    MatrixXd eKva(as<Map<MatrixXd> >(Kva));
    VectorXd ey(as<Map<MatrixXd> >(y));
    MatrixXd eX(as<Map<MatrixXd> >(X));

    VectorXd result = getMLsoln(hsq, eKva, ey, eX, reml);

    int p = result.size()-3;
    NumericVector beta(p);
    for(int i=0; i<p; i++)
        beta[i] = result[i+3];

    return List::create(Named("sigsq") =     result[0],
                        Named("rss") =       result[1],
                        Named("logdetXSX") = result[2],
                        Named("beta") =      beta);
}

// calcLL
// calculate log likelihood for fixed value of hsq
// sigmasq = total variance = sig^2_g + sig^2_e
//
// hsq   = heritability
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = boolean indicating whether to use REML (vs ML)
// logdetXpX = log det X'X; if NA, it's calculated
VectorXd calcLL(double hsq, VectorXd Kva, VectorXd y,
                MatrixXd X, bool reml=true, double logdetXpX=NA_REAL)
{
    int n = Kva.size();
    int p = X.cols();

    // estimate beta and sigma^2
    VectorXd ml_soln = getMLsoln(hsq, Kva, y, X, reml);
    double sigsq = ml_soln[0];
    double rss   = ml_soln[1];
    double logdetXSX = ml_soln[2];

    // calculate log likelihood
    double loglik = (double)n*log(rss);
    for(int i=0; i<n; i++)
        loglik += log(hsq*Kva[i] + 1.0 - hsq);
    loglik *= -0.5;

    if(reml) {
        if(NumericVector::is_na(logdetXpX)) { // need to calculate it
            MatrixXd XpX(calc_xpx(X));
            std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);
            logdetXpX=0.0;
            for(int i=0; i<p; i++) logdetXpX += log(e.first[i]);
        }

        loglik += 0.5*(p*log(2*M_PI*sigsq) + logdetXpX - logdetXSX);
    }

    VectorXd result(p+2);
    result[0] = loglik;
    result[1] = sigsq;
    for(int i=0; i<p; i++) result[i+2] = ml_soln[i+3];
    return result;
}

// calcLL (version called from R)
// [[Rcpp::export]]
List R_calcLL(double hsq, NumericVector Kva, NumericVector y,
              NumericMatrix X, bool reml=true, double logdetXpX=NA_REAL)
{
    MatrixXd eKva(as<Map<MatrixXd> >(Kva));
    VectorXd ey(as<Map<MatrixXd> >(y));
    MatrixXd eX(as<Map<MatrixXd> >(X));

    VectorXd result = calcLL(hsq, eKva, ey, eX, reml, logdetXpX);

    int p = result.size()-2;
    NumericVector beta(p);
    for(int i=0; i<p; i++)
        beta[i] = result[i+2];

    return List::create(Named("loglik") =    result[0],
                        Named("sigsq") =     result[1],
                        Named("beta") =      beta);
}
