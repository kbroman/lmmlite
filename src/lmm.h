// linear mixed model via RcppEigen
#ifndef LMM_H
#define LMM_H

MatrixXd calc_XpX(const MatrixXd& X);

// eigen decomp
List eigen_decomp(const NumericMatrix &A);

#endif // LMM_H
