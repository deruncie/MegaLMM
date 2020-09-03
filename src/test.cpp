// #include <math.h>
// #include <iostream>
// #include "MegaLMM_types.h"
// 
// using namespace Eigen;
// 
// MatrixXf rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec(n*p);
//   for(int i = 0; i < n*p; i++){
//     X_vec[i] = ziggr.norm();
//   }
//   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
//   return(X_mat.cast<float>());
// }
// struct General_Matrix_f2 {
//   MatrixXf dense;
//   SpMat sparse;
//   bool triangular;
//   bool isDense;
//   bool isNULL = false;
//   General_Matrix_f2(MatrixXf dense_, SpMat sparse_,bool triangular_, bool isDense_,bool isNULL_) : dense(dense_), sparse(sparse_), triangular(triangular_), isDense(isDense_), isNULL(isNULL_){}
//   MatrixXf solve(MatrixXf) const;
//   // General_Matrix_f solve(General_Matrix_f) const;
//   MatrixXf tsolve(MatrixXf) const;
//   MatrixXf operator*(MatrixXf) const;
//   MatrixXf crossprod(MatrixXf) const;
//   float get_log_det();
//   int rows() const;
//   int cols() const;
// };
// // Loads a sparse or dense matrix passed from R into a \code{General_Matrix_f2} object
// General_Matrix_f2 load_General_Matrix_f2(SEXP X_, bool triangular) {
//   MatrixXf null_d = MatrixXf::Zero(0,0);
//   if(Rf_isNull(X_)) {
//     SpMat null_s = null_d.sparseView();
//     General_Matrix_f2 Xm(null_d,null_s,triangular,false,true);
//     return(Xm);
//   } else if(Rf_isMatrix(X_)){
//     MatrixXf X = as<MatrixXf >(X_);
//     SpMat null_s = null_d.sparseView();
//     General_Matrix_f2 Xm(X,null_s,triangular,true,false);
//     return(Xm);
//   } else{
//     SpMat X = as<SpMat>(X_);
//     General_Matrix_f2 Xm(null_d,X,triangular,false,false);
//     return(Xm);
//   }
// }
// 
// MatrixXf General_Matrix_f2::solve(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense.triangularView<Upper>().solve(Y));
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse.triangularView<Upper>().solve(Y));
//     }
//   } else{
//     // Note: these Cholesky's could be stored so they could be re-used.
//     if(isDense) {
//       return(dense.ldlt().solve(Y));
//     } else{
//       Eigen::SimplicialLDLT<SpMat> ldlt(sparse);
//       return(ldlt.solve(Y));
//     }
//   }
// }
// 
// // MatrixXf General_Matrix_f2::solve(General_Matrix_f2 Y) const {
// //   if(isNULL) return(Y);
// //   if(triangular) {
// //     if(isDense) {
// //       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
// //       if(Y.isDense) {
// //         return(dense.triangularView<Upper>().solve(Y.dense));
// //       } else {
// //         MatrixXf Z(dense.rows,Y.cols);
// //         for(int j = 0; j < Y.cols; j++) Z.col(j) = dense.triangularView<Upper>().solve(Y.sparse.col(j));
// //       }
// //     } else{
// //       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
// //       if(Y.isDense) {
// //         return(sparse.triangularView<Upper>().solve(Y.dense));
// //       } else {
// //         return(sparse.triangularView<Upper>().solve(Y.sparse));
// //       }
// //     }
// //   } else{
// //     // Note: these Cholesky's could be stored so they could be re-used.
// //     if(isDense) {
// //       if(Y.isDense) {
// //         return(dense.ldlt().solve(Y.dense));
// //       } else {
// //         return(dense.ldlt().solve(Y.sparse));
// //       }
// //     } else{
// //       Eigen::SimplicialLDLT<SpMat> ldlt(sparse);
// //       if(Y.isDense) {
// //         return(ldlt.solve(Y.dense));
// //       } else {
// //         return(ldlt.solve(Y.sparse));
// //       }
// //     }
// //   }
// // }
// MatrixXf General_Matrix_f2::tsolve(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense.transpose().triangularView<Lower>().solve(Y));
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse.transpose().triangularView<Lower>().solve(Y));
//     }
//   } else{
//     // Note: these Cholesky's could be stored so they could be re-used.
//     if(isDense) {
//       return(dense.transpose().ldlt().solve(Y));
//     } else{
//       Eigen::SimplicialLDLT<SpMat> ldlt(sparse.transpose());
//       return(ldlt.solve(Y));
//     }
//   }
// }
// 
// MatrixXf General_Matrix_f2::operator*(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense.triangularView<Upper>() * Y);
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse * Y);
//     }
//   } else {
//     if(isDense) {
//       if(Y.rows() != dense.cols()) stop("Wrong dimension for Y");
//       return(dense * Y);
//     } else{
//       if(Y.rows() != sparse.cols()) stop("Wrong dimension for Y");
//       return(sparse * Y);
//     }
//   }
// }
// MatrixXf General_Matrix_f2::crossprod(MatrixXf Y) const {
//   if(isNULL) return(Y);
//   if(triangular) {
//     if(isDense) {
//       if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
//       return(dense.transpose().triangularView<Lower>() * Y);
//     } else{
//       if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
//       return(sparse.transpose().triangularView<Lower>() * Y);
//     }
//   } else{
//     if(isDense) {
//       if(Y.rows() != dense.rows()) stop("Wrong dimension for Y");
//       return(dense.transpose() * Y);
//     } else{
//       if(Y.rows() != sparse.rows()) stop("Wrong dimension for Y");
//       return(sparse.transpose() * Y);
//     }
//   }
// }
// float General_Matrix_f2::get_log_det() {
//   if(isNULL) stop("no determinant of null matrix");
//   if(!triangular) stop("not implemented for non-triangular matrix");
//   if(isDense) {
//     return(dense.diagonal().array().log().sum());
//   } else {
//     float log_det = 0;
//     for(int j = 0; j < sparse.rows(); j++) {
//       log_det += std::log(sparse.coeffRef(j,j));
//     }
//     return(log_det);
//   }
// }
// int General_Matrix_f2::rows() const {
//   int r;
//   if(isDense) {
//     r = dense.rows(); 
//   } else {
//     r = sparse.rows();
//   }
//   return(r);
// }
// int General_Matrix_f2::cols() const {
//   int c;
//   if(isDense) {
//     c = dense.cols(); 
//   } else {
//     c = sparse.cols();
//   }
//   return(c);
// }
// 
// // Code to convert list of R matrices (sparse or dense) into a thread-safe object
// //
// // @param X_list List of matrices (each can be dgCMatrix or matrix)
// // @param X_vector \code{std::vector} of \code{General_Matrix_f2} type which will be populated
// void load_General_Matrix_f2_list(const Rcpp::List X_list, std::vector<General_Matrix_f2>& X_vector, bool triangular){
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     X_vector.push_back(load_General_Matrix_f2(Xi_, triangular));
//   }
// }
// 
// 
// 
// // [[Rcpp::export()]]
// VectorXf tsf(SpMat chol_ZtZ_Kinv, VectorXf b) {
//   return(chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b));
// }
// 
// // [[Rcpp::export()]]
// VectorXd tsd(SpMatd chol_ZtZ_Kinv, VectorXd b) {
//   return(chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b));
// }
//   
// 
// VectorXf sample_MME_single_diagR2(
//     VectorXf y,           // nx1
//     General_Matrix_f2 Z,    // nxr dgCMatrix or dense
//     SpMat chol_ZtZ_Kinv,       // rxr CsparseMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
//     float tot_Eta_prec,   // float
//     float pe,            // float
//     VectorXf randn_theta  // rx1
// ){
//   VectorXf b;
//   b = Z.crossprod(y)*pe;
//   Rcout << b.sum() << " ";
//   // if(Rf_isMatrix(Z_)) {
//   //   MatrixXf Z = as<MatrixXf>(Z_);
//   //   b = Z.transpose() * y * pe;
//   // } else{
//   //   SpMat Z = as<SpMat>(Z_);
//   //   b = Z.transpose() * y * pe;
//   // }
//   // if(Z.isDense) {
//   //   b = Z.dense.transpose() * y * pe;
//   // } else{
//   //   b = Z.sparse.transpose() * y * pe;
//   // }
//   b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b);
//   b = b / sqrt(tot_Eta_prec);
//   Rcout << b.sum() << " ";
//   b += randn_theta;
//   Rcout << b.sum() << " ";
//   b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b / sqrt(tot_Eta_prec));
//   Rcout << b.sum() << "\n";
//   return(b);
// }
// 
// 
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// // E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// // For complete data, ie no missing obs.
// // [[Rcpp::export()]]
// MatrixXf sample_MME_ZKZts_c2(
//     MatrixXf Y,                    // nxp
//     SEXP Z_,
//     VectorXf tot_Eta_prec,         // px1
//     Rcpp::List chol_ZtZ_Kinv_list_,      // List or R st RtR = ZtZ_Kinv
//     MatrixXf h2s,                  // n_RE x p
//     VectorXi h2s_index                 // px1
// ) {
//   
//   General_Matrix_f2 Z = load_General_Matrix_f2(Z_,false);
//   
//   int p = Y.cols();
//   int r = Z.cols();
//   // if(Z.isDense) {
//   //   r = Z.dense.cols();
//   // } else{
//   //   r = Z.sparse.cols();
//   // }
//   
//   MatrixXf randn_theta = rstdnorm_mat2(r,p);
//   
//   std::vector<General_Matrix_f2> chol_ZtZ_Kinv_list;
//   load_General_Matrix_f2_list(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list, true);
//   
//   MatrixXf U(r,p);
//   ArrayXf h2_e = 1.0 - h2s.colwise().sum().array();
//   ArrayXf pes = tot_Eta_prec.array() / h2_e.array();
//   
// // #pragma omp parallel for
//   for(std::size_t j = 0; j < p; j++){
//     Rcout << j << " ";
//     int h2_index = h2s_index[j] - 1;
//     // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
//     U.col(j) = sample_MME_single_diagR2(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index].sparse, tot_Eta_prec[j], pes[j],randn_theta.col(j));
//   }
//   
//   return(U);
// }
// 
// 
// // [[Rcpp::export()]]
// List LDLt2(SEXP A_) {
//   if(Rf_isMatrix(A_)){
//     MatrixXd A = as<MatrixXd >(A_);
//     Eigen::LDLT<MatrixXd> ldlt_A;
//     ldlt_A.compute(A);
//     MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
//     MatrixXd P = ldlt_A.transpositionsP() * I;
//     VectorXd d = ldlt_A.vectorD();
//     MatrixXd L = ldlt_A.matrixL();
//     SpMatd Lsp = L.sparseView();
//     if(static_cast<float>(Lsp.nonZeros()) / I.size() > 0.25) {
//       return(List::create(
//           Named("P") = P.sparseView(),
//           Named("L") = L,
//           Named("d") = d
//       ));
//     } else{
//       return(List::create(
//           Named("P") = P.sparseView(),
//           Named("L") = Lsp,
//           Named("d") = d
//       ));
//     }
//   } else{
//     SpMatd A = as<SpMatd>(A_);
//     Eigen::SimplicialLDLT<SpMatd> ldlt_A;
//     ldlt_A.compute(A);
//     MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
//     MatrixXd P = ldlt_A.permutationP() * I;
//     return(List::create(
//         Named("P") = P.sparseView(),
//         Named("L") = ldlt_A.matrixL(),
//         Named("d") = ldlt_A.vectorD()
//     ));
//   }
// }
