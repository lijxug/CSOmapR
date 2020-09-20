#define EIGEN_USE_BLAS
#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
// #include <RcppArmadillo.h>
#include <RcppEigen.h>

//#include <RcppParallel.h>

using Eigen::MatrixXd;
using Eigen::VectorXd;
// using arma::mat;
using namespace Rcpp;


// [[Rcpp::export]]
MatrixXd calc_d_rcpp(MatrixXd ydata){
  MatrixXd y_square = ydata.array().square();
  VectorXd sum_ydata = y_square.rowwise().sum();
  
  MatrixXd YYt = ydata * ydata.transpose();
  
  MatrixXd d = -2 * YYt;
  d.rowwise() += sum_ydata.transpose();
  d.colwise() += sum_ydata;
  
  return(d);
}

// [[Rcpp::export]]
MatrixXd update_grads_rcpp(MatrixXd grads, MatrixXd ydata, MatrixXd stiffnesses){
  for(int i = 0; i<grads.rows(); i++){
    MatrixXd ydata_add = -1 * ydata;
    //VectorXd ydata_chosenrow = ydata.block(i, 0, 1, ydata.cols()).transpose();
    VectorXd ydata_chosenrow = ydata.row(i).transpose();
    ydata_add.rowwise() += ydata_chosenrow.transpose();
    
    VectorXd stiff_chosencol = stiffnesses.col(i);
    ydata_add.array().colwise() *= stiff_chosencol.array();
    
    VectorXd y_colsum = ydata_add.colwise().sum();
    grads.row(i) = y_colsum.transpose();
  }
  return(grads);
}


