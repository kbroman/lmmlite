// linear mixed model via RcppEigen
#ifndef LMM_H
#define LMM_H

struct lmm_fit {
    double hsq;
    VectorXd beta;
    double sigmasq;
    double loglik;
    double rss;
    double logdetXSX;
};

struct calcLL_args {
    VectorXd Kva;
    VectorXd y;
    MatrixXd X;
    bool reml;
    double logdetXpX;
};


// calc X'X
MatrixXd calc_XpX(const MatrixXd& X);

// calc X'X (version to be called from R)
NumericMatrix R_calc_xpx(const NumericMatrix& X);

// eigen decomposition
//    returns eigenvalues and transposed eigenvectors
std::pair<Eigen::VectorXd, Eigen::MatrixXd> eigen_decomp(Eigen::MatrixXd A);

// eigen decomposition
//    returns list with eigenvalues and transposed eigenvectors
List R_eigen_decomp(const NumericMatrix &A);

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
                         MatrixXd X, bool reml);

// getMLsoln (version called from R)
List R_getMLsoln(double hsq, NumericVector Kva, NumericVector y,
                 NumericMatrix X, bool reml);

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
                      MatrixXd X, bool reml, double logdetXpX);

// calcLL (version called from R)
List R_calcLL(double hsq, NumericVector Kva, NumericVector y,
              NumericMatrix X, bool reml, double logdetXpX);

// just the negative log likelihood, for the optimization
double negLL(double x, struct calcLL_args *args);

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
                      bool reml, bool check_boundary,
                      double logdetXpX, double tol);

// fitLMM (version called from R)
List R_fitLMM(NumericVector Kva, NumericVector y, NumericMatrix X,
              bool reml, bool check_boundary,
              double logdetXpX, double tol);

#endif // LMM_H
