#include <math.h>
#include <iostream>
#include "MegaLMM_types.h"

using namespace Eigen;


RMatrix load_RMatrix(SEXP X_, bool triangular) {
  MatrixXd null_d = MatrixXd::Zero(0,0);
  if(Rf_isNull(X_)) {
    SpMat null_s = null_d.sparseView();
    MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
    Map<MatrixXd> M_null_d(null_d.data(),0,0);
    RMatrix Xm(M_null_d,M_null_s,triangular,false,true);
    return(Xm);
  } else if(Rf_isMatrix(X_)){
    Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
    SpMat null_s = null_d.sparseView();
    MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
    RMatrix Xm(X,M_null_s,triangular,true,false);
    return(Xm);
  } else{
    MSpMat X = as<MSpMat>(X_);
    Map<MatrixXd> M_null_d(null_d.data(),0,0);
    RMatrix Xm(M_null_d,X,triangular,false,false);
    return(Xm);
  }
}

MatrixXd RMatrix::solve(MatrixXd Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense.triangularView<Upper>().solve(Y));
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse.triangularView<Upper>().solve(Y));
    }
  } else{
    // Note: these Cholesky's could be stored so they could be re-used.
    if(isDense) {
      return(dense.ldlt().solve(Y));
    } else{
      Eigen::SimplicialLDLT<SpMat> ldlt(sparse);
      return(ldlt.solve(Y));
    }
  }
}

MatrixXd RMatrix::tsolve(MatrixXd Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense.transpose().triangularView<Lower>().solve(Y));
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse.transpose().triangularView<Lower>().solve(Y));
    }
  } else{
    // Note: these Cholesky's could be stored so they could be re-used.
    if(isDense) {
      return(dense.transpose().ldlt().solve(Y));
    } else{
      Eigen::SimplicialLDLT<SpMat> ldlt(sparse.transpose());
      return(ldlt.solve(Y));
    }
  }
}

MatrixXd RMatrix::operator*(MatrixXd Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense.triangularView<Upper>() * Y);
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse * Y);
    }
  } else {
    if(isDense) {
      if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
      return(dense * Y);
    } else{
      if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
      return(sparse * Y);
    }
  }
}
MatrixXd RMatrix::crossprod(MatrixXd Y) const {
  if(isNULL) return(Y);
  if(triangular) {
    if(isDense) {
      if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
      return(dense.transpose().triangularView<Lower>() * Y);
    } else{
      if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
      return(sparse.transpose().triangularView<Lower>() * Y);
    }
  } else{
    if(isDense) {
      if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
      return(dense.transpose() * Y);
    } else{
      if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
      return(sparse.transpose() * Y);
    }
  }
}
float RMatrix::get_log_det() {
  if(isNULL) stop("no determinant of null matrix");
  if(!triangular) stop("not implemented for non-triangular matrix");
  if(isDense) {
    return(dense.diagonal().array().log().sum());
  } else {
    float log_det = 0;
    for(int j = 0; j < sparse.rows(); j++) {
      log_det += std::log(sparse.coeffRef(j,j));
    }
    return(log_det);
  }
}
int RMatrix::rows() const {
  int r;
  if(isDense) {
    r = dense.rows(); 
  } else {
    r = sparse.rows();
  }
  return(r);
}
int RMatrix::cols() const {
  int c;
  if(isDense) {
    c = dense.cols(); 
  } else {
    c = sparse.cols();
  }
  return(c);
}


