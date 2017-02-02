#include <RcppEigen.h>
#include "sigmoid.hpp"
#include <limits>

typedef std::numeric_limits<double> double_lim;
//[[Rcpp::export]]
Eigen::ArrayXd betavar(const Eigen::ArrayXd &p,const Eigen::ArrayXd &mu,Eigen::ArrayXd &s){
  return p*(s+(1-p)*mu.square());
}
//[[Rcpp::export]]
double intklbeta_rssbvsr(const Eigen::ArrayXd &alpha,const Eigen::ArrayXd &mu,const Eigen::ArrayXd &sigma_square,const double &sigma_beta_square){
  double tres = alpha.sum()+alpha.matrix().transpose()*((sigma_square/sigma_beta_square).log()).matrix();
  double sres = alpha.matrix().transpose()*(sigma_square+mu.square()).matrix();
  double thres = alpha.matrix().transpose()*((alpha+double_lim::epsilon()).log()).matrix();
  double fres = (1-alpha).matrix().transpose()*(1-alpha+double_lim::epsilon()).log().matrix();
  return (tres-sres/sigma_beta_square)*0.5-thres-fres;
}

//[[Rcpp::export]]
double intgamma(const Eigen::ArrayXd &logodds, const Eigen::ArrayXd &alpha){
  Eigen::ArrayXd tres = ((alpha-1)*logodds+logsigmoid(logodds));
  return tres.sum();
}