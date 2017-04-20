#include <RcppEigen.h>
#include "rssr.h"
//#include "mkl.h"
#include <math.h>

//[[Rcpp::export(name="logsigmoid")]]
Eigen::ArrayXd logsigmoid_exp(const arrayxd_external x){
  return(logsigmoid(x));
  //  return -(-x).unaryExpr(&logpexp);
}


//[[Rcpp::export(name="sigmoid")]]
Eigen::ArrayXd sigmoid_exp(const arrayxd_external x){
  return(sigmoid(x));
  //  return -(-x).unaryExpr(&logpexp);
}



