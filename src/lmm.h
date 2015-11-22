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


#endif // LMM_H
