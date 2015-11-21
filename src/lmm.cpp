// linear mixed model via RcppEigen

// [[Rcpp::depends(RcppEigen)]]

#include <RcppEigen.h>

using namespace Rcpp;
using namespace Eigen;

#include "lmm.h"

// calc X'X
Eigen::MatrixXd calc_XpX(const Eigen::MatrixXd& X)
{
    int n = X.cols();

    return MatrixXd(n,n).setZero().selfadjointView<Lower>()
        .rankUpdate(X.transpose());
}

// eigen decomp
// [[Rcpp::export]]
List eigen_decomp(NumericMatrix A)
{
    MatrixXd AA(as<Map<MatrixXd> >(A));
    const Eigen::SelfAdjointEigenSolver<MatrixXd> VLV(AA);
    return List::create(Named("vectors") = VLV.eigenvectors(),
                        Named("values") = VLV.eigenvalues());
}
