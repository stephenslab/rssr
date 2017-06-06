#include <RcppEigen.h>
#include "rssr.h"
#include <cstdio>
#include <cmath>
#include <RcppParallel.h>
#include <algorithm>
#include <progress.hpp>
// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]


using namespace tbb;

Eigen::ArrayXd initialize_s_m(const c_arrayxd_internal sesquare, const double sigma_beta_square){
  return((sesquare*(sigma_beta_square))/(sesquare+(sigma_beta_square)));
}



Eigen::ArrayXd initialize_ssrat_m(const c_arrayxd_internal s, const double sigma_beta_square){
  return((s/sigma_beta_square).log());
}


 double calculate_mtp_m(const c_arrayxd_internal alpha_r,const c_arrayxd_internal alpha_v,const c_arrayxd_internal mu_r, const c_arrayxd_internal mu_v){
  
  return(  -std::sqrt(((alpha_r).square()).sum()+((mu_r).square()).sum())/std::sqrt(((alpha_v).square()).sum()+((mu_v).square()).sum()));
}


class rssr{

private:
  const size_t p;
  const size_t n;


  const bool use_squarem;
  const bool use_lnztol;
  const int itermax;
  const double tolerance;

  const Eigen::ArrayXd logodds;
  const Eigen::ArrayXd sigma_beta;

  const size_t logodds_size;
  const size_t sigb_size;
  const size_t gridsize;
  const ArrayXd sigma_beta_square;
  
  const arrayxd_external betahat;
  const Eigen::ArrayXd sesquare;
  
  const Eigen::ArrayXd  q;
  Eigen::ArrayXXd s;
  Eigen::ArrayXXd ssrat;
  
  
  const Matrix_external SiRiS;
  
    
  std::vector<double> lnZ_d;
  std::vector<double> lnZ0_d;
  std::vector<double> lnZ00_d;
  
  
  
  std::vector<double> alpha_d;
  std::vector<double> mu_d;
  std::vector<double> SiRiSr_d;
  
  std::vector<double> alpha0_d;
  std::vector<double> mu0_d;
  std::vector<double> SiRiSr0_d;
  
  std::vector<double> alpha1_d;
  std::vector<double> mu1_d;
  std::vector<double> SiRiSr1_d;

  std::vector<double> mtp_d;

  
  
  std::vector<double> alpha_r;
  std::vector<double> alpha_v;
  
  std::vector<double> mu_r;
  std::vector<double> mu_v;

  bool reverse;
  std::vector<double>* final_max_err;
  std::vector<double>* pvevec;
  Eigen::ArrayXi final_iternum;
  std::vector<int> forward_range;
  std::vector<int> reverse_range;

  const size_t max_backtrack=12;
  
public:
  rssr(const Matrix_external _SiRiS,
       const arrayxd_external _sigma_beta,
       const arrayxd_external _logodds,
       const arrayxd_external _betahat,
       const arrayxd_external _se,
       const arrayxd_external _alpha,
       const arrayxd_external _mu,
       const arrayxd_external _SiRiSr,
       const double _tolerance,
       const int _itermax,
       const bool _lnz_tol,
       const bool _use_squarem,
       const size_t _n) : p(_betahat.size()),
			  n(_n),
			  use_squarem(_use_squarem),
			  use_lnztol(_lnz_tol),
			  itermax(_itermax),
			  tolerance(_tolerance),
			  logodds(_logodds),
			  sigma_beta(_sigma_beta),
			  logodds_size(logodds.size()),
			  sigb_size(sigma_beta.size()),
			  gridsize(sigb_size),
			  sigma_beta_square(sigma_beta*sigma_beta),
			  betahat(_betahat),
			  sesquare(_se*_se),
			  q(betahat/sesquare),
			  SiRiS(_SiRiS),
			  alpha(p,gridsize),
			  mu(p,gridsize),
			  SiRiSr(p,gridsize),
			  final_iternum(p),
			  final_max_err(p),
			  pvevec(p),
			  alpha0(alpha),
			  mu0(mu),
			  SiRiSr1(SiRiSr),
			  alpha1(alpha),
			  mu1(mu),
			  alpha_r(alpha),
			  alpha_v(alpha),
			  mu_r(mu),
			  mu_v(mu),
			  lnZ(gridsize),
			  lnZ0(lnZ),
			  lnZ00(lnZ0)
  {

    

    mtp=-1;
    for(int i=0;i<p;i++){
      forward_range[i]=i;
      reverse_range[i]=p-1-i;
      final_max_err(i)=1;
      //      iternum[i]=0;
      }
    for(const int& i: forward_range){
      alpha.col(i)=_alpha;
      mu.col(i)=_mu;
      SiRiSr.col(i)=_SiRiSr;
    }
    for(size_t j=0;j<gridsize;j++){
      s.col(j)=initialize_s_m(sesquare,sigma_beta_square(j));
      ssrat.col(j)=initialize_ssrat_m(s.col(j),sigma_beta_square(j));
    }
    
  }
  void operator()(const blocked_range<size_t> &r)const{

    double max_err=1;
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;
   

    
    calc_lnZ(r);
    int iternum=0;
    
    
      
    lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
    while(max_err>tolerance){
      lnZ00.segment(rstart,rsize)=lnZ0.segment(rstart,rsize);
      lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
      alpha0.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
      mu0.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);

      bool reverse = iternum%2!=0;
      rss_vb_iter(r,reverse);
      alpha1.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
      mu1.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);

      rss_vb_iter(r,reverse);

      alpha_r.block(0,rstart,p,rsize)=(alpha1.block(0,rstart,p,rsize)-alpha0.block(0,rstart,p,rsize));
      alpha_v.block(0,rstart,p,rsize)=(alpha.block(0,rstart,p,rsize)-alpha1.block(0,rstart,p,rsize)-alpha_r.block(0,rstart,p,rsize));

      mu_r.block(0,rstart,p,rsize)=(mu1.block(0,rstart,p,rsize)-mu0.block(0,rstart,p,rsize));
      mu_v.block(0,rstart,p,rsize)=(mu.block(0,rstart,p,rsize)-mu1.block(0,rstart,p,rsize)-mu_r.block(0,rstart,p,rsize));

      squarem_step(r);

      rss_vb_iter(r,reverse);

      calc_lnZ(r);

      squarem_backtrack(r);
      calc_err(r,max_err);
      if(iternum>itermax){
	for(size_t  j=r.begin();j!=r.end();j++){
	  final_max_err(j)=max_err;
	  final_iternum(j)=iternum;
	}
	compute_pve(r);
	break;
      }
      iternum=iternum+1;
    }
    for(size_t  j=r.begin();j!=r.end();j++){
      final_max_err(j)=max_err;
      final_iternum(j)=iternum;
    }
    compute_pve(r);
  }
    
    void calc_lnZ(const blocked_range<size_t> & r)const {
      
      for(size_t j=r.begin(); j!=r.end();j++){
	lnZ(j)=calculate_lnZ(q,alpha.col(j)*mu.col(j),SiRiSr.col(j),
			     logodds(j),sesquare,
			     alpha.col(j),mu.col(j),s.col(j),sigma_beta(j));
      }
    }

    void rss_vb_iter(const blocked_range<size_t> & r,bool reverse)const{
    
    

      for(size_t j=r.begin(); j!=r.end();j++){
	//        bool reverse = iternum[j]%2!=0;				  
	if(reverse){
	  for(const int& i: reverse_range){
	    rss_varbvsr_update(betahat.coeffRef(i),
			       sesquare.coeffRef(i),
			       sigma_beta_square(j),
			       SiRiS.col(i),
			       s(i,j),
			       SiRiSr.col(j),
			       SiRiSr(i,j),
			       ssrat(i,j),
			       logodds(j),
			       alpha(i,j),
			       mu(i,j));
	  }
	}else{
	  for(const int& i: forward_range){
	    rss_varbvsr_update(betahat.coeffRef(i),
			       sesquare.coeffRef(i),
			       sigma_beta_square(j),
			       SiRiS.col(i),
			       s(i,j),
			       SiRiSr.col(j),
			       SiRiSr(i,j),
			       ssrat(i,j),
			       logodds(j),
			       alpha(i,j),
			       mu(i,j));
	  }
	}
      }
    }

    void squarem_step(const blocked_range<size_t> & r)const {
      
      for(size_t j=r.begin(); j!=r.end();j++){
	double tmtp =calculate_mtp_m(alpha_r,alpha_v,mu_r,mu_v);
	mtp(j)= tmtp;
	  if(tmtp>(-1)){
	    tmtp=-1;
	    mtp(j) =-1;
	  }else{
	    if(!std::isnan(tmtp)){
	      alpha.col(j)=alpha0.col(j)-2*tmtp*alpha_r.col(j)+(tmtp*tmtp)*(alpha_v.col(j));
	      mu.col(j)=mu0.col(j)-2*tmtp*(mu_r.col(j))+(tmtp*tmtp)*(mu_v.col(j));
	      SiRiSr.col(j)=(SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	    }
	  }
      }
    }
      
    void squarem_backtrack(const blocked_range<size_t> & r)const {
      
      for(size_t j=r.begin();j!=r.end();j++){
	double tmtp=mtp[j];
	double tlnZ=lnZ[j];
	double tlnZ0=lnZ0[j];
	  if(tmtp<(-1) && (tlnZ < tlnZ0)){
	    size_t num_bt=0;
	    while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	      tmtp=0.5*(tmtp-1);
	      alpha.col(j) = alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*(alpha.col(j)-2*alpha1.col(j)+alpha0.col(j));
	      mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	      SiRiSr = (SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	      rss_vb_iter(blocked_range<size_t>(j,j),false);
	      rss_vb_iter(blocked_range<size_t>(j,j),true);
	      
	      calc_lnZ(blocked_range<size_t>(j,j));
	      num_bt=num_bt+1;
	    }
	    if(num_bt==max_backtrack){
	      alpha.col(j)=alpha0.col(j);
	      mu.col(j)=mu0.col(j);
	      SiRiSr.col(j) = (SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	      calc_lnZ(blocked_range<size_t>(j,j));
	      lnZ0(j)=lnZ00(j);
	      // Rcpp::Rcerr<<"Final!: num_bt:"<<num_bt<<" step_size:"<<mtp<<"lnZ:"<<lnZ<<" lnZ0:"<<lnZ0<<std::endl;
	    }
	  }
	mtp[j]=tmtp;
      }
    }
    
    void calc_err(const blocked_range<size_t> &r,double &max_err)const {
      
      max_err=0;
      if(use_lnztol){
	
	for(size_t j=r.begin(); j!=r.end();j++){
	  double tmax_err=rel_err(lnZ(j),lnZ0(j));
	  if(tmax_err>max_err){
	    max_err=tmax_err;
	  }
	}
      }else{
	for(size_t j=r.begin(); j!=r.end();j++){
	  double tmax_err=find_maxerr(alpha.col(j),alpha0.col(j),alpha.col(j)*mu.col(j),alpha0.col(j)*mu0.col(j));
	  if(tmax_err>max_err){
	    max_err=tmax_err;
	  }
	}
      }
    }
  void compute_pve(const blocked_range<size_t> &r)const {
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;
    //    for(size_t j=r.begin();j!=r.end();j++){
    pvevec.segment(rstart,rsize)=SiRiSr.block(0,rstart,p,rsize).matrix().transpose()*(alpha.block(0,rstart,p,rsize)*mu.block(0,rstart,p,rsize)).matrix();
  }
  
  Rcpp::DataFrame give_results(){
    using namespace Rcpp;
    return(Rcpp::DataFrame::create(_["logodds"]=Rcpp::wrap(logodds),
				   _["sigb"]=Rcpp::wrap(sigma_beta),
				   _["rel_err"]=Rcpp::wrap(final_max_err),
				   _["iterations"]=Rcpp::wrap(final_iternum),
				   _["alpha_mean"]=Rcpp::wrap(alpha.colwise().mean()),
				   _["pve"]=Rcpp::wrap(pvevec),
				   _["lnZ"]=Rcpp::wrap(lnZ)));
  }
    
};



//[[Rcpp::export]]
Rcpp::DataFrame grid_search_rss_varbvsr_alt(
    const  Matrix_external SiRiS,
    const arrayxd_external sigma_beta,
    const arrayxd_external logodds,
    const arrayxd_external  betahat,
    const arrayxd_external  se,
    const arrayxd_external talpha0,
    const arrayxd_external tmu0,
    const arrayxd_external tSiRiSr0,
    const double tolerance,
    const int itermax,
    const int n,
    Rcpp::LogicalVector lnz_tol){
  //  std::cout<<"Starting grid_search_rss_varbvsr"<<std::endl;

  bool islnz_tol = lnz_tol(0);

  if(sigma_beta.minCoeff()<=0){
    Rcpp::stop("sigma_beta must be strictly positive");
  }
  size_t sigb_size= sigma_beta.size();
  size_t logodds_size=logodds.size();
  size_t tot_size=sigb_size;
  if(tot_size!=logodds_size){
    Rcpp::stop("Length of sigma_beta must equal length of logodds");
  }
  rssr rssr_obj(SiRiS,sigma_beta,logodds,betahat,se,talpha0,tmu0,tSiRiSr0,tolerance,itermax,islnz_tol,true,n);
  parallel_for(blocked_range<size_t>(0,tot_size),rssr_obj);
  
  return(rssr_obj.give_results());
} 
