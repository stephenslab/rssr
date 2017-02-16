#include <RcppEigen.h>
#include "sigmoid.hpp"
#include "kl.hpp"


//[[Rcpp::export]]
Eigen::ArrayXd betavar(const Eigen::ArrayXd &p,const Eigen::ArrayXd &mu,const Eigen::ArrayXd &s){
  return p*(s+(1-p)*mu.square());
}
//[[Rcpp::export]]
double intklbeta_rssbvsr(const Eigen::ArrayXd &alpha,const Eigen::ArrayXd &mu,const Eigen::ArrayXd &sigma_square, double sigma_beta_square){
  double tres = alpha.sum()+alpha.matrix().transpose()*((sigma_square/sigma_beta_square).log()).matrix();
  double sres = alpha.matrix().transpose()*(sigma_square+mu.square()).matrix();
  double thres = alpha.matrix().transpose()*((alpha+double_lim::epsilon()).log()).matrix();
  double fres = (1-alpha).matrix().transpose()*(1-alpha+double_lim::epsilon()).log().matrix();
  return (tres-sres/sigma_beta_square)*0.5-thres-fres;
}

//[[Rcpp::export]]
double intgamma(double logodds, const Eigen::ArrayXd &alpha){
  Eigen::ArrayXd tres = ((alpha-1)*logodds+logsigmoid(logodds));
  return tres.sum();
}

//[[Rcpp::export]]
double rel_err(double p0,double p1){
  return fabs(p0-p1)/(fabs(p0)+fabs(p1)+double_lim::epsilon());
}

//[[Rcpp::export]]
double find_maxerr(const Eigen::ArrayXd &alpha,
                   const Eigen::ArrayXd &alpha0,
                   const Eigen::ArrayXd &r,
                   const Eigen::ArrayXd &r0){
  double terr=0;
  double cerr;
  size_t p=alpha.size();
  for(size_t i=0;i<p;i++){
    if(fabs(alpha(i))>1e-6){
      cerr=fabs(alpha(i)-alpha0(i))/(fabs(alpha(i))+fabs(alpha0(i))+double_lim::epsilon());
      if(cerr>terr){
        terr=cerr;
      }
    }
    if(fabs(r(i))>1e-6){
      cerr=fabs(r(i)-r0(i))/(fabs(r(i))+fabs(r0(i))+double_lim::epsilon());
    }
    if(cerr>terr){
      terr=cerr;
    }
  }
  return(terr);
}


//[[Rcpp::export]]
double calculate_lnZ(const Eigen::VectorXd &q,
                     const Eigen::VectorXd &r,
                     const Eigen::VectorXd &SiRiSr,
                     double logodds,
                     const Eigen::VectorXd &sesquare,
                     const Eigen::VectorXd &alpha,
                     const Eigen::VectorXd &mu,
                     const Eigen::VectorXd &s,
                     double sigb){
  

  double lnz0 = q.dot(r)-0.5*r.dot(SiRiSr)+intgamma(logodds,alpha.array());
  double lnz1 = lnz0-0.5*(1/sesquare.array()).matrix().dot(betavar(alpha.array(),mu.array(),s.array()).matrix());
  double lnz2 = lnz1+intklbeta_rssbvsr(alpha.array(),mu.array(),s.array(),sigb*sigb);
  return(lnz2);
}







