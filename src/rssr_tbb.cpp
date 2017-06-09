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
using namespace RcppParallel;

Eigen::ArrayXd initialize_s_m(const c_arrayxd_internal sesquare, const double sigma_beta_square){
  return((sesquare*(sigma_beta_square))/(sesquare+(sigma_beta_square)));
}



Eigen::ArrayXd initialize_ssrat_m(const c_arrayxd_internal s, const double sigma_beta_square){
  return((s/sigma_beta_square).log());
}


 double calculate_mtp_m(const c_arrayxd_internal alpha_r,const c_arrayxd_internal alpha_v,const c_arrayxd_internal mu_r, const c_arrayxd_internal mu_v){
  
  return(  -std::sqrt(((alpha_r).square()).sum()+((mu_r).square()).sum())/std::sqrt(((alpha_v).square()).sum()+((mu_v).square()).sum()));
}


struct rssr: public Worker{


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
  
    
  Eigen::ArrayXd lnZ;
  Eigen::ArrayXd lnZ0;
  Eigen::ArrayXd lnZ00;
  
  
  
  Eigen::ArrayXXd alpha;
  Eigen::ArrayXXd mu;
  Eigen::ArrayXXd SiRiSr;
  
  Eigen::ArrayXd alpha0;
  Eigen::ArrayXd mu0;
  Eigen::ArrayXd SiRiSr0;
  
  Eigen::ArrayXd alpha1;
  Eigen::ArrayXd mu1;
  Eigen::ArrayXd SiRiSr1;

  Eigen::ArrayXd mtp;

  
  
  Eigen::ArrayXXd alpha_r;
  Eigen::ArrayXXd alpha_v;
  
  Eigen::ArrayXXd mu_r;
  Eigen::ArrayXXd mu_v;

  Eigen::ArrayXd rvec;

  Eigen::ArrayXd final_max_err;
  Eigen::ArrayXd pvevec;
  Eigen::ArrayXd alpha_mean;
  
  
  Eigen::ArrayXi final_iternum;
  std::vector<int> forward_range;
  std::vector<int> reverse_range;

  const size_t max_backtrack=12;
  
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
			  s(p,gridsize),
			  ssrat(p,gridsize),
			  SiRiS(_SiRiS),
			  alpha(p,gridsize),
			  mu(p,gridsize),
			  SiRiSr(p,gridsize),
			  final_iternum(gridsize),
			  final_max_err(gridsize),
			  pvevec(gridsize),
			  alpha0(alpha),
			  mu0(mu),
			  SiRiSr1(SiRiSr),
			  alpha1(alpha),
			  mu1(mu),
			  alpha_r(alpha),
			  alpha_v(alpha),
			  mu_r(mu),
			  mu_v(mu),
			  mtp(gridsize),
			  alpha_mean(gridsize),
			  lnZ(gridsize),
			  lnZ0(lnZ),
			  rvec(gridsize),
			  lnZ00(lnZ0)
  {

    
    
    forward_range.resize(p);
    reverse_range.resize(p);
    final_max_err.resize(gridsize);
    //    mtp=-1;
    // Rcpp::Rcout<<"p is :"<<p<<std::endl;
    // Rcpp::Rcout<<"gridsize is :"<<gridsize<<std::endl;
    
    for(int i=0;i<p;i++){
      forward_range[i]=i;
      reverse_range[i]=i;

      //      iternum[i]=0;
      }
    std::reverse(std::begin(reverse_range),std::end(reverse_range));
    // Rcpp::Rcout<<"Starting second loop! :"<<forward_range.size()<<std::endl;
    for(size_t j=0;j<gridsize;j++){
      alpha.col(j)=_alpha;
      mu.col(j)=_mu;
      SiRiSr.col(j)=_SiRiSr;
      rvec(j)=0;
      final_max_err(j)=1;
    }
    
    // //<<"Starting third loop! :"<<gridsize<<std::endl;
    for(size_t j=0;j<gridsize;j++){
      s.col(j)=initialize_s_m(sesquare,sigma_beta_square(j));
      ssrat.col(j)=initialize_ssrat_m(s.col(j),sigma_beta_square(j));
    }



    

    if(p!=alpha.rows()){
      Rcpp::stop("Shape of alpha is wrong");
    }
    if(p!=mu.rows()){
      Rcpp::stop("Shape of mu is wrong");
    }
    if(p!=mu.rows()){
      Rcpp::stop("Shape of mu is wrong");
    }
    // //<<"Finished third loop!"<<std::endl;
    
  }
  void operator()(std::size_t begin, std::size_t end){
    
    blocked_range<size_t> r(begin,end);
    double max_err=1;
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;
   

    
    calc_lnZ(r);
    int iternum=0;
    //<<"Initializing lnZ Segment "<<std::endl;    
    calc_lnZ(r);
    //<<"Initializing lnZ0 Segment "<<std::endl;
    lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
    while(max_err>tolerance){
      if(iternum==0)
      //<<"Setting lnZ00 Segment "<<std::endl;
      lnZ00.segment(rstart,rsize)=lnZ0.segment(rstart,rsize);
      //<<"Setting lnZ00 Segment "<<std::endl;
      lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
      //<<"Setting alpha0 block "<<std::endl;
      alpha0.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
      //<<"Setting mu0 block "<<std::endl;
      mu0.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);

      bool reverse = iternum%2!=0;
      //<<"Performing first Iteration "<<std::endl;
      rss_vb_iter(r,reverse);
      //<<"Setting alpha1 block "<<std::endl;
      alpha1.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
      mu1.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
      //<<"Performing second Iteration "<<std::endl;
      rss_vb_iter(r,reverse);

      alpha_r.block(0,rstart,p,rsize)=(alpha1.block(0,rstart,p,rsize)-alpha0.block(0,rstart,p,rsize));
      alpha_v.block(0,rstart,p,rsize)=(alpha.block(0,rstart,p,rsize)-alpha1.block(0,rstart,p,rsize)-alpha_r.block(0,rstart,p,rsize));

      mu_r.block(0,rstart,p,rsize)=(mu1.block(0,rstart,p,rsize)-mu0.block(0,rstart,p,rsize));
      mu_v.block(0,rstart,p,rsize)=(mu.block(0,rstart,p,rsize)-mu1.block(0,rstart,p,rsize)-mu_r.block(0,rstart,p,rsize));
      //<<"Squarem Step "<<std::endl;
      squarem_step(r);
      //<<"Third Iteration "<<std::endl;
      rss_vb_iter(r,reverse);

      calc_lnZ(r);
      //<<"Squarem Backtrack "<<std::endl;
      squarem_backtrack(r);
      //<<"Error Calculation "<<std::endl;
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
    //<<"Error Calculation "<<std::endl;
    for(size_t  j=r.begin();j!=r.end();j++){
      final_max_err(j)=max_err;
      final_iternum(j)=iternum;
      alpha_mean(j)=alpha.col(j).mean();
    }
    compute_pve(r);

  }
    
    void calc_lnZ(const blocked_range<size_t> & r) {
      
      for(size_t j=r.begin(); j!=r.end();j++){
	lnZ(j)=calculate_lnZ(q,alpha.col(j)*mu.col(j),SiRiSr.col(j),
			     logodds(j),sesquare,
			     alpha.col(j),mu.col(j),s.col(j),sigma_beta(j));
      }
    }

    void rss_vb_iter(const blocked_range<size_t> & r,bool reverse){
    
    size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;
    
    

      // for(size_t j=r.begin(); j!=r.end();j++){
	//        bool reverse = iternum[j]%2!=0;				  
	if(reverse){
	  for(const int& i: reverse_range){
	      
	      
	      rvec.segment(rstart,rsize) = alpha.row(i).segment(rstart,rsize) * mu.row(i).segment(rstart,rsize);
	    

	    mu.row(i).segment(rstart,rsize) = s.row(i).segment(rstart,rsize) * (q(i) + rvec.segment(rstart,rsize)/sesquare(i) - SiRiSr.row(i).segment(rstart,rsize));
	    alpha.row(i).segment(rstart,rsize) = (1+(-(logodds.segment(rstart,rsize) + 0.5 * (ssrat.row(i).segment(rstart,rsize) + mu.row(i).segment(rstart,rsize).square()/s.row(i).segment(rstart,rsize)))).exp()).inverse();
	    
	    
	    for(size_t j=r.begin();j!=r.end();j++){
	      SiRiSr.col(j)+=SiRiS.col(i).array()*(alpha(i,j)*mu(i,j)-rvec(j));
	    }
//	    SiRiSr.block(0,rstart,p,rsize).array().rowwise() += ((alpha.row(i).segment(rstart,rsize)*mu.row(i).segment(rstart,rsize)-rvec.segment(rstart,rsize)).colwise()*SiRiS.col(i).array()).array();
	    
	    }

	}else{
	  for(const int& i: forward_range){
	    rvec.segment(rstart,rsize) = alpha.row(i).segment(rstart,rsize) * mu.row(i).segment(rstart,rsize);
	    
	    // mu.row(i).segment(rstart,rsize);
	    mu.row(i).segment(rstart,rsize) = s.row(i).segment(rstart,rsize) * (q(i) + rvec.segment(rstart,rsize)/sesquare(i) - SiRiSr.row(i).segment(rstart,rsize));
	    alpha.row(i).segment(rstart,rsize) = (1+(-(logodds.segment(rstart,rsize) + 0.5 * (ssrat.row(i).segment(rstart,rsize) + mu.row(i).segment(rstart,rsize).square()/s.row(i).segment(rstart,rsize)))).exp()).inverse();
	    
	    for(size_t j=r.begin();j!=r.end();j++){
	      SiRiSr.col(j)+=SiRiS.col(i).array()*(alpha(i,j)*mu(i,j)-rvec(j));
	    }
//	    SiRiSr.block(0,rstart,p,rsize).array().rowwise() += ((alpha.row(i).segment(rstart,rsize)*mu.row(i).segment(rstart,rsize)-rvec.segment(rstart,rsize)).colwise()*SiRiS.col(i).array()).array();
	    
	  }
	}
    }
  
  void squarem_step(const blocked_range<size_t> & r) {
      
      for(size_t j=r.begin(); j!=r.end();j++){
	
	//<<"mtp_calc:";
	double tmtp =calculate_mtp_m(alpha_r.col(j),alpha_v.col(j),mu_r.col(j),mu_v.col(j));
	//<<tmtp<<std::endl;
	mtp(j)= tmtp;
	  if(tmtp>(-1)){
	    tmtp=-1;
	    mtp(j) =-1;
	  }else{
	    if(!std::isnan(tmtp)){
	      //<<"alpha_step:"<<std::endl;
	      alpha.col(j)=alpha0.col(j)-2*tmtp*alpha_r.col(j)+(tmtp*tmtp)*(alpha_v.col(j));
	      //<<"mu_step:"<<std::endl;
	      mu.col(j)=mu0.col(j)-2*tmtp*(mu_r.col(j))+(tmtp*tmtp)*(mu_v.col(j));
	      //<<"SiRiSr_step:"<<std::endl;
	      SiRiSr.col(j)=(SiRiS*(alpha.col(j)*mu.col(j)).matrix()).array();
	    }
	  }
      }
    }
      
    void squarem_backtrack(const blocked_range<size_t> & r) {
      
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
    
    void calc_err(const blocked_range<size_t> &r,double &max_err) {
      
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
  void compute_pve(const blocked_range<size_t> &r) {
    // size_t rstart=r.begin();
    // size_t rsize= r.end()-rstart;
    for(size_t j=r.begin();j!=r.end();j++){
      pvevec(j)=(SiRiSr.col(j).matrix().transpose()*(alpha.col(j)*mu.col(j)).matrix());
      pvevec(j)/=n;
    }
  }
  
  Rcpp::DataFrame give_results(){
    using namespace Rcpp;
    Rcpp::NumericVector rlogodds=Rcpp::wrap(logodds);
    Rcpp::NumericVector rsigb=Rcpp::wrap(sigma_beta);
    Rcpp::NumericVector rrel_err=Rcpp::wrap(final_max_err);
    Rcpp::IntegerVector riterations=Rcpp::wrap(final_iternum);

    Rcpp::NumericVector ralphamean=Rcpp::wrap(alpha_mean);
    Rcpp::NumericVector rpve=Rcpp::wrap(pvevec);
    Rcpp::NumericVector rlnz=Rcpp::wrap(lnZ);

    std::vector<size_t> sizes(7);
    sizes.push_back(rlogodds.size());
    sizes.push_back(rsigb.size());
	
    sizes.push_back(rrel_err.size());
    sizes.push_back(riterations.size());
    sizes.push_back(ralphamean.size());
    sizes.push_back(rpve.size());
    sizes.push_back(rlnz.size());
    if ( std::adjacent_find( sizes.begin(), sizes.end(), std::not_equal_to<int>() ) == sizes.end() ){
      Rcpp::stop("Return vectors aren't all the same size");
    }
    
    return(Rcpp::DataFrame::create(_["logodds"]=rlogodds,
				   _["sigb"]=rsigb,
				   _["rel_err"]=rrel_err,
				   _["iterations"]=riterations,
				   _["alpha_mean"]=ralphamean,
				   _["pve"]=rpve,
				   _["lnZ"]=rlnz));
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
    Rcpp::LogicalVector lnz_tol,const int grainsize=1){
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
  Rcpp::Rcout  <<"Initializing rssr_obj"<<std::endl;
  rssr rssr_obj(SiRiS,sigma_beta,logodds,betahat,se,talpha0,tmu0,tSiRiSr0,tolerance,itermax,islnz_tol,true,n);
  Rcpp::Rcout  <<"Starting Serial Algorithm!"<<std::endl;
  rssr_obj(0,tot_size);
  //    parallelFor(0,tot_size,rssr_obj,grainsize);
  Rcpp::Rcout<<"Finished!"<<std::endl;
  return(rssr_obj.give_results());
} 









// struct rssr_alt: public Worker{
//   
//   Eigen::MatrixXd 

// 
// Rcpp::DataFrame grid_rss_varbvsr_alt2(
//     const c_Matrix_internal SiRiS,
//     const c_arrayxd_internal sigma_beta,
//     const c_arrayxd_internal logodds,
//     const c_arrayxd_internal  betahat,
//     const c_arrayxd_internal  se,
//     const c_arrayxd_internal talpha0,
//     const c_arrayxd_internal tmu0,
//     const c_arrayxd_internal tSiRiSr0,
//     double tolerance,
//     int itermax,
//     bool isVerbose,
//     bool islnz_tol){
//   
//   //  std::cout<<"Starting grid_rss_varbvsr (tbb)"<<std::endl;
//   
//   using namespace Rcpp;
//   size_t sigb_size= sigma_beta.size();
//   size_t logodds_size=logodds.size();
//   
//   size_t tot_size=sigb_size;
//   if(tot_size!=logodds_size){
//     Rcpp::stop("Length of sigma_beta must equal length of logodds");
//   }
//   
//   Eigen::ArrayXd npivec(tot_size);
//   Eigen::ArrayXd nlzvec(tot_size);
//   Eigen::ArrayXd errvec(tot_size);
//   Eigen::ArrayXi itervec(tot_size);
//   Eigen::ArrayXd pvevec(tot_size);
//   
//   
//   Rcpp::LogicalVector verbose(1);
//   verbose(0)=isVerbose;
//   Rcpp::LogicalVector lnz_tol(1);
//   lnz_tol(0)=islnz_tol;
//   // static affinity_partitioner ap;
//   parallel_for(blocked_range<size_t>(0,tot_size),
//                [&](const blocked_range<size_t>& r)  {
//                  for(size_t t=r.begin(); t!=r.end(); t++){
//                    size_t i=t;
//                    size_t j=t;
//                    //                   sigbvec(t)=sigma_beta(j);
//                    Eigen::ArrayXd copy_alpha(talpha0);
//                    Eigen::ArrayXd copy_mu(tmu0);
//                    Eigen::ArrayXd copy_SiRiSr(tSiRiSr0);
//                    //                   lovec(t)=logodds(i);
//                    int iter=0;
//                    double max_err=1;
//                    double retvec=rss_varbvsr_squarem_iter(SiRiS,
//                                                           sigma_beta(j),
//                                                           logodds(i),
//                                                           betahat,
//                                                           se,
//                                                           copy_alpha,
//                                                           copy_mu,
//                                                           copy_SiRiSr,
//                                                           tolerance,
//                                                           itermax,
//                                                           lnz_tol,iter,max_err);
//                    nlzvec[t]=retvec;
//                    itervec[t]=iter;
//                    errvec[t]=max_err;
//                    npivec[t]=copy_alpha.mean();
//                    pvevec[t]=copy_SiRiSr.matrix().transpose()*(copy_alpha*copy_mu).matrix();
//                  }});
//   
//   //  Rcpp::Rcout<<"mean lnZ is: "<<mean(nlzvec)<<std::endl;
//   return(Rcpp::DataFrame::create(_["logodds"]=logodds,
//                                  _["sigb"]=sigma_beta,
//                                  _["rel_err"]=errvec,
//                                  _["iterations"]=itervec,
//                                  _["alpha_mean"]=npivec,
//                                  _["pve"]=pvevec,
//                                  _["lnZ"]=nlzvec));
// }





