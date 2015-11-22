// linear mixed model via RcppEigen
#ifndef LMM_H
#define LMM_H

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
// reml  = boolean indicating whether to use REML (or ML)
VectorXd getMLsoln(double hsq, VectorXd Kva, VectorXd y,
                   MatrixXd X, bool reml);

// getMLsoln (version called from R)
List R_getMLsoln(double hsq, NumericVector Kva, NumericVector y,
                 NumericMatrix X, bool reml);

#endif // LMM_H
