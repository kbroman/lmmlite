// linear mixed model via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

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
