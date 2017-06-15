#include <RcppEigen.h>
#include "rssr.h"



Eigen::ArrayXd betavar(const c_arrayxd_internal p,const c_arrayxd_internal mu,const c_arrayxd_internal s){
  return p*(s+(1-p)*mu.square());
}

//[[Rcpp::export(name="betavar")]]
Eigen::ArrayXd exp_betavar(const arrayxd_external p, const arrayxd_external mu, const arrayxd_external s){
  return betavar(p,mu,s);
}



double intklbeta_rssbvsr(const c_arrayxd_internal alpha,const c_arrayxd_internal mu,const c_arrayxd_internal sigma_square, double sigma_beta_square){
  double tres = alpha.sum()+alpha.matrix().transpose()*((sigma_square/sigma_beta_square).log()).matrix();
  double sres = alpha.matrix().transpose()*(sigma_square+mu.square()).matrix();
  double thres = alpha.matrix().transpose()*((alpha+double_lim::epsilon()).log()).matrix();
  double fres = (1-alpha).matrix().transpose()*(1-alpha+double_lim::epsilon()).log().matrix();
  return (tres-sres/sigma_beta_square)*0.5-thres-fres;
}

//[[Rcpp::export(name="intklbeta_rssbvsr")]]
double exp_intklbeta_rssbvsr(const arrayxd_external alpha,const arrayxd_external mu,const arrayxd_external sigma_square, double sigma_beta_square){
  return(intklbeta_rssbvsr(alpha,mu,sigma_square,sigma_beta_square));
}





// This is equivalent to sum(alpha*log(pi)+(1-alpha)*log(1-pi))
double intgamma(double logodds, const c_arrayxd_internal alpha){
  Eigen::ArrayXd tres = ((alpha-1)*logodds+logsigmoid(logodds));
  return tres.sum();
}

//[[Rcpp::export(name="intgamma")]]
double exp_intgamma(double logodds, const arrayxd_external alpha){
  return(intgamma(logodds,alpha));
}




//[[Rcpp::export]]
double rel_err(double p0,double p1){
  return fabs(p0-p1)/(fabs(p0)+fabs(p1)+double_lim::epsilon());
}

Eigen::ArrayXd rel_err(c_arrayxd_internal p0,c_arrayxd_internal p1){
  return ((p0-p1).cwiseAbs())/(p0.cwiseAbs()+p1.cwiseAbs()+double_lim::epsilon());
}



double find_maxerr(const c_arrayxd_internal alpha,
                   const c_arrayxd_internal alpha0,
                   const c_arrayxd_internal r,
                   const c_arrayxd_internal r0){
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

//[[Rcpp::export(name="find_maxerr")]]
double exp_find_maxerr(const arrayxd_external alpha,
                   const arrayxd_external alpha0,
                   const arrayxd_external r,
                   const arrayxd_external r0){
  return(find_maxerr(alpha,alpha0,r,r0));
}


double update_logodds(const c_arrayxd_internal alpha){
  double pi=alpha.mean();
  return log((pi+double_lim::epsilon())/((1-pi)+double_lim::epsilon()));
}

//[[Rcpp::export(name="update_logodds")]]
double exp_update_logodds(const arrayxd_external alpha){
  return(update_logodds(alpha));
}
  


// double calculate_lnZ(const double* qp,
// 		     const double* SiRiSrp,
// 		     const double logodds,
// 		     const double* sesquarep,
// 		     const double* alphap,
// 		     const double* mup,
// 		     const double* sp,
// 		     const double sigb,const size_t p){
  
  
//   mdvec q(qp,p);
//   mdvec alpha(alphap,p);
//   mdvec mu(alphap,p);
//   double lnz0 = q.dot(r)-0.5*r.dot(SiRiSr)+intgamma(logodds,alpha.array());
//   // if(!std::isfinite(lnz0)){
//   //   Rcpp::stop("lnZ0 is not finite");
//   // }
//   double lnz1 = lnz0-0.5*(1/sesquare.array()).matrix().dot(betavar(alpha.array(),mu.array(),s.array()).matrix());
//   // if(!std::isfinite(lnz0)){
//   //   Rcpp::stop("lnZ1 is not finite");
//   // }
//   double lnz2 = lnz1+intklbeta_rssbvsr(alpha.array(),mu.array(),s.array(),sigb*sigb);
  
  
//   // if(!std::isfinite(lnz0)){
//   //   Rcpp::stop("lnZ2 is not finite");
//   // }
//   return(lnz2);
// }



  double calculate_lnZ(const c_vectorxd_internal q,
                       const c_vectorxd_internal r,
                       const c_vectorxd_internal SiRiSr,
                       const double logodds,
                       const c_vectorxd_internal sesquare,
                       const c_vectorxd_internal alpha,
                       const c_vectorxd_internal mu,
                       const c_vectorxd_internal s,
                       const double sigb){
    
    
    double lnz0 = q.dot(r)-0.5*r.dot(SiRiSr)+intgamma(logodds,alpha.array());
    // if(!std::isfinite(lnz0)){
    //   Rcpp::stop("lnZ0 is not finite");
    // }
    double lnz1 = lnz0-0.5*(1/sesquare.array()).matrix().dot(betavar(alpha.array(),mu.array(),s.array()).matrix());
    // if(!std::isfinite(lnz0)){
    //   Rcpp::stop("lnZ1 is not finite");
    // }
    double lnz2 = lnz1+intklbeta_rssbvsr(alpha.array(),mu.array(),s.array(),sigb*sigb);
    // if(!std::isfinite(lnz0)){
    //   Rcpp::stop("lnZ2 is not finite");
    // }
    return(lnz2);
  }





// double exp_calculate_lnZ(const vectorxd_external q,
//                      const vectorxd_external r,
//                      const vectorxd_external SiRiSr,
//                      double logodds,
//                      const vectorxd_external sesquare,
//                      const vectorxd_external alpha,
//                      const vectorxd_external mu,
//                      const vectorxd_external s,
//                      double sigb){
// 
// 
//   return(calculate_lnZ(q,r,SiRiSr,logodds,sesquare,alpha,mu,s,sigb));
// }



//[[Rcpp::export(name="calculate_lnZ")]]
double wrap_calculate_lnZ(const arrayxd_external alpha_mu,
                                             const Matrix_external SiRiS,
                                             const double sigma_beta,
                                             const double logodds,
                                             const arrayxd_external betahat,
                                             const arrayxd_external se){
  
  
  size_t p=betahat.size();
  // const arrayxd_external alpha,
  // const arrayxd_external mu,
  // const arrayxd_external SiRiSr,
  using namespace Eigen;  
  
  Eigen::ArrayXd talpha=alpha_mu.head(p);
  Eigen::ArrayXd tmu=alpha_mu.tail(p);
  Eigen::ArrayXd tSiRiSr=SiRiS*(talpha*tmu).matrix();
//  Eigen::ArrayXd ret_alpha_mu(alpha_mu.size());
  double sigma_beta_square=sigma_beta*sigma_beta;
  Eigen::ArrayXd sesquare=se.square();
  Eigen::ArrayXd  s= (sesquare*(sigma_beta*sigma_beta))/(sesquare+(sigma_beta*sigma_beta));
  Eigen::ArrayXd ssrat((s/sigma_beta_square).log());
  Eigen::ArrayXd  q= betahat/sesquare;
  return(calculate_lnZ(q,talpha*tmu,tSiRiSr,logodds,sesquare,talpha,tmu,s,sigma_beta));
  
}


