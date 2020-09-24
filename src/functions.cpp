#define EIGEN_USE_BLAS
// [[Rcpp::depends(RcppEigen)]]
#include <Rcpp.h>
#include <RcppEigen.h>


// using arma::mat;
using namespace Rcpp;
using Eigen::MatrixXd;
using Eigen::VectorXd;


//' Calculate d using Rcpp
//' @param ydata A matrix of coordinates
//' @param v2 Second value
//' @return d
// [[Rcpp::export]]
Eigen::MatrixXd calc_d_rcpp(Eigen::MatrixXd ydata){
  Eigen::MatrixXd y_square = ydata.array().square();
  Eigen::VectorXd sum_ydata = y_square.rowwise().sum();
  
  Eigen::MatrixXd YYt = ydata * ydata.transpose();
  
  Eigen::MatrixXd d = -2 * YYt;
  d.rowwise() += sum_ydata.transpose();
  d.colwise() += sum_ydata;
  
  return(d);
}

//' Update gradients using Rcpp
//' @param grads A matrix of gradients
//' @param ydata A matrix of coordinates
//' @param stiffnesses A matrix of stiffnesses
//' @return updated grads
// [[Rcpp::export]]
Eigen::MatrixXd update_grads_rcpp(Eigen::MatrixXd grads, Eigen::MatrixXd ydata, Eigen::MatrixXd stiffnesses){
  for(int i = 0; i<grads.rows(); i++){
    Eigen::MatrixXd ydata_add = -1 * ydata;
    //Eigen::VectorXd ydata_chosenrow = ydata.block(i, 0, 1, ydata.cols()).transpose();
    Eigen::VectorXd ydata_chosenrow = ydata.row(i).transpose();
    ydata_add.rowwise() += ydata_chosenrow.transpose();
    
    Eigen::VectorXd stiff_chosencol = stiffnesses.col(i);
    ydata_add.array().colwise() *= stiff_chosencol.array();
    
    Eigen::VectorXd y_colsum = ydata_add.colwise().sum();
    grads.row(i) = y_colsum.transpose();
  }
  return(grads);
}


