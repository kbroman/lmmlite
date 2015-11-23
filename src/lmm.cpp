// linear mixed model via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include <math.h>
#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "brent_fmin.h"
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


// eigen + rotation
// perform eigen decomposition of kinship matrix
// and rotate phenotype and covariate matrices by transpose of eigenvectors
struct eigenrot eigen_rotation(MatrixXd K, MatrixXd y, MatrixXd X)
{
    std::pair<VectorXd,MatrixXd> e = eigen_decomp(K);
    MatrixXd yrot = e.second * y;
    MatrixXd Xrot = e.second * X;

    struct eigenrot result;
    result.Kva = e.first;
    result.Kve = e.second;
    result.y = yrot;
    result.X = Xrot;

    return result;
}

// eigen + rotation
// [[Rcpp::export]]
List R_eigen_rotation(NumericMatrix K, NumericMatrix y, NumericMatrix X)
{
    MatrixXd KK(as<Map<MatrixXd> >(K));
    MatrixXd yy(as<Map<MatrixXd> >(y));
    MatrixXd XX(as<Map<MatrixXd> >(X));

    struct eigenrot result = eigen_rotation(KK, yy, XX);

    return List::create(Named("Kva") = result.Kva,
                        Named("Kve") = result.Kve,
                        Named("y") = result.y,
                        Named("X") = result.X);
}

// calculate log det X'X
double calc_logdetXpX(MatrixXd X)
{
    MatrixXd XpX(calc_xpx(X)); // calc X'X
    int p = X.cols();

    // eigen decomposition of X'X
    std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);

    // calculate log det X'X
    double result=0.0;
    for(int i=0; i<p; i++) result += log(e.first[i]);

    return result;
}

// calculate log det X'X (version to be called from R)
// [[Rcpp::export]]
double R_calc_logdetXpX(NumericMatrix X)
{
    MatrixXd XX(as<Map <MatrixXd> >(X));

    return calc_logdetXpX(XX);
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
struct lmm_fit getMLsoln(double hsq, VectorXd Kva, VectorXd y,
                   MatrixXd X, bool reml=true)
{
    const int n = Kva.size();
    const int p = X.cols();
    struct lmm_fit result;

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
    result.rss = rss(0,0);
    result.sigmasq = result.rss/(double)(n-p);
    result.beta = beta.col(0);
    result.logdetXSX = logdetXSX; // determinant (if REML)

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

    struct lmm_fit result = getMLsoln(hsq, eKva, ey, eX, reml);

    return List::create(Named("sigmasq") =     result.sigmasq,
                        Named("rss") =       result.rss,
                        Named("logdetXSX") = result.logdetXSX,
                        Named("beta") =      result.beta);
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
struct lmm_fit calcLL(double hsq, VectorXd Kva, VectorXd y,
                MatrixXd X, bool reml=true, double logdetXpX=NA_REAL)
{
    int n = Kva.size();
    int p = X.cols();

    // estimate beta and sigma^2
    struct lmm_fit ml_soln = getMLsoln(hsq, Kva, y, X, reml);

    // calculate log likelihood
    double loglik = (double)n*log(ml_soln.rss);
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

        loglik += 0.5*(p*log(2 * M_PI * ml_soln.sigmasq) + logdetXpX - ml_soln.logdetXSX);
    }

    ml_soln.loglik = loglik;
    return ml_soln;
}

// calcLL (version called from R)
// [[Rcpp::export]]
List R_calcLL(double hsq, NumericVector Kva, NumericVector y,
              NumericMatrix X, bool reml=true, double logdetXpX=NA_REAL)
{
    MatrixXd eKva(as<Map<MatrixXd> >(Kva));
    VectorXd ey(as<Map<MatrixXd> >(y));
    MatrixXd eX(as<Map<MatrixXd> >(X));

    struct lmm_fit result = calcLL(hsq, eKva, ey, eX, reml, logdetXpX);

    return List::create(Named("loglik") =    result.loglik,
                        Named("sigmasq") =   result.sigmasq,
                        Named("beta") =      result.beta);
}

// just the negative log likelihood, for the optimization
double negLL(double x, struct calcLL_args *args)
{
    struct lmm_fit result = calcLL(x, args->Kva, args->y, args->X,
                                   args->reml, args->logdetXpX);

    return -result.loglik;
}

// fitLMM
// Optimize log liklihood over hsq
//
// Kva   = eigenvalues of kinship matrix
// y     = rotated vector of phenotypes
// X     = rotated matrix of covariates
// reml  = boolean indicating whether to use REML (vs ML)
// check_boundary = if true, explicity check 0.0 and 1.0 boundaries
// logdetXpX = log det X'X; if NA, it's calculated
// tol   = tolerance for convergence
struct lmm_fit fitLMM(VectorXd Kva, VectorXd y, MatrixXd X,
                      bool reml=true, bool check_boundary=true,
                      double logdetXpX=NA_REAL, double tol=1e-4)
{
    struct lmm_fit result;

    // calculate log det XpX, if necessary
    // (note same befor and after it's "rotated" by eigenvec of kinship matrix
    if(reml && NumericVector::is_na(logdetXpX)) {
        MatrixXd XpX(calc_xpx(X));
        std::pair<VectorXd, MatrixXd> e = eigen_decomp(XpX);
        int p = X.cols();
        logdetXpX=0.0;
        for(int i=0; i<p; i++) logdetXpX += log(e.first[i]);
    }

    // function arguments for calcLL
    struct calcLL_args args;
    args.Kva = Kva;
    args.y = y;
    args.X = X;
    args.reml = reml;
    args.logdetXpX = logdetXpX;

    double hsq = qtl2_Brent_fmin(0.0, 1.0, (double (*)(double, void*)) negLL, &args, tol);
    result = calcLL(hsq, Kva, y, X, reml, logdetXpX);
    result.hsq = hsq;

    if(check_boundary) {
        struct lmm_fit boundary_result;
        boundary_result = calcLL(0.0, Kva, y, X, reml, logdetXpX);
        if(boundary_result.loglik > result.loglik) {
            result = boundary_result;
            result.hsq = 0.0;
        }
        boundary_result = calcLL(1.0, Kva, y, X, reml, logdetXpX);
        if(boundary_result.loglik > result.loglik) {
            result = boundary_result;
            result.hsq = 1.0;
        }
    }

    return result;
}

// fitLMM (version called from R)
// [[Rcpp::export]]
List R_fitLMM(NumericVector Kva, NumericVector y, NumericMatrix X,
              bool reml=true, bool check_boundary=true,
              double logdetXpX=NA_REAL, double tol=1e-4)
{
    MatrixXd eKva(as<Map<MatrixXd> >(Kva));
    VectorXd ey(as<Map<MatrixXd> >(y));
    MatrixXd eX(as<Map<MatrixXd> >(X));

    struct lmm_fit result = fitLMM(eKva, ey, eX, reml, check_boundary,
                                   logdetXpX, tol);

    return List::create(Named("loglik") =    result.loglik,
                        Named("hsq") =       result.hsq,
                        Named("sigmasq") =   result.sigmasq,
                        Named("beta") =      result.beta);
}
