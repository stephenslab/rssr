#include <RcppEigen.h>
#include "rssr.h"

// [[Rcpp::depends(RcppProgress)]]
//[[Rcpp::depends(RcppParallel)]]





# define rassert(EX) (void)((EX) || (__assert (#EX, __FILE__, __LINE__),0))
void __assert (const char *msg, const char *file, int line) {
  char buffer [200];
  snprintf( buffer, 200, "Assert: %s at %s line #%d\n", msg, file,
	    line );
  ::Rf_error( buffer );
}

template<typename T, typename D> void zero_out(T &dyn_array, const D fillval){
#ifndef EIGEN_NO_DEBUG
  std::fill(dyn_array.local().begin(),dyn_array.local().end(),fillval);
#endif
}
  
  

using namespace Eigen;


Fit_wrap::Fit_wrap(double* fit_,double* fit0_,double* fit1_,const size_t p_, const size_t gridsize_): p(p_),
												      gridsize(gridsize_),
												      tls(false),
												      fit_p(fit_),
												      fit0_p(fit0_),
												      fit1_p(fit1_),
												      fitl(tbbdvec(0)),
												      fitl0(tbbdvec(0)),
												      fitl1(tbbdvec(0))
  
  {
#ifndef NDEBUG
    m2darray fit(fit_p,p,gridsize);
    m2darray fit0(fit0_p,p,gridsize);
    m2darray fit1(fit1_p,p,gridsize);
    rassert(fit.cols()==gridsize &&"fit should have gridsize no of columns");
    rassert(fit0.cols()==gridsize &&"fit0 should have gridsize no of columns");
    rassert(fit1.cols()==gridsize &&"fit1 should have gridsize no of columns");
    rassert(isnan(fit).sum()==0 && "fit has no NaN  values");
#endif
    
  }


Fit_wrap::Fit_wrap(double* fit_,const size_t p_, const size_t gridsize_): p(p_),
									  gridsize(gridsize_),
									  fit_p(fit_),
									  tls(false),
									  fit0_p(NULL),
									  fit1_p(NULL),
									  fitl(tbbdvec(0)),
									  fitl0(tbbdvec(0)),
									  fitl1(tbbdvec(0)){

  
  
}



Fit_wrap::Fit_wrap(double* fit_,double* fit0_, double* fit1_,const size_t p_):p(p_),
									      gridsize(1),
									      fit_p(fit_),
									      tls(true),
									      fit0_p(fit0_),
									      fit1_p(fit1_),
									      fitl(tbbdvec(0)),
									      fitl0(tbbdvec(0)),
									      fitl1(tbbdvec(0)){


  
}

Fit_wrap::Fit_wrap(double* fit_,const size_t p_): p(p_),
						  gridsize(1),
						  fit_p(fit_),
						  tls(true),
						  fit0_p(NULL),
						  fit1_p(NULL),
						  fitl(tbbdvec(0)),
						  fitl0(tbbdvec(0)),
						  fitl1(tbbdvec(0)){

  
  
}


void Fit_wrap::tlresize(const size_t rsize){
  

  fitl.local().reserve(rsize*p);
  mdarray ofit(fit_p,p);
  m2darray nfit(fitl.local().data(),p,rsize);
  for(int i=0; i<rsize;i++){
    nfit.col(i)=ofit;
  }
  rassert(nfit.cols()==rsize &&"fit should have same num of cols as rsize");
  rassert(isnan(nfit).sum()==0 && "fit has no NaN  values");
//  fit_p=fitl.local().data();
  
  
  if(fit0_p!=NULL){
    //    Rcpp::Rcout<<"resizing fit0"<<std::endl;
    mdarray ofit0(fit0_p,p);
    fitl0.local().resize(rsize*p);
    m2darray nfit0(fitl0.local().data(),p,rsize);
    
    for(int i=0; i<rsize;i++){
      nfit0.col(i)=ofit;
    }
//    fit0_p=fitl0.local().data();
    rassert(nfit0.cols()==rsize &&"fit0 should have rsize no of columns");
    rassert(isnan(nfit0).sum()==0 && "fit0 has no NaN  values");
  }
  else{
    //    // Rcpp::Rcout<<"NOT resizing fit0"<<std::endl;
  }
  if(fit1_p!=NULL){
    // Rcpp::Rcout<<"resizing fit1"<<std::endl;
    mdarray ofit1(fit1_p,p);
    fitl1.local().resize(rsize*p);
    m2darray nfit1(fitl1.local().data(),p,rsize);
    
    for(int i=0; i<rsize;i++){
      nfit1.col(i)=ofit1;
    }
    rassert(nfit1.cols()==rsize &&"fit1 should have rsize no of columns");
    rassert(isnan(nfit1).sum()==0 && "fit1 has no NaN  values");
//    fit1_p=fitl1.local().data();
  }
  else{
    // Rcpp::Rcout<<"NOT resizing fit1"<<std::endl;
  }
}
		   
Siris::Siris(double* siris_,const size_t p_):siris_p(siris_),p(p_){
    
    //rassert(SiRiS.cols()==SiRiS.rows() && "SiRiS should be square");
    //rassert(SiRiS.cols()==p && "SiRiS should have p cols");
    //rassert(isnan(SiRiS.array()).sum()==0 && "sesquare has no NaN  values");

}
Siris::Siris(Eigen::MatrixXd &siris_):siris_p(siris_.data()),p(siris_.cols()){

  //rassert(SiRiS.cols()==SiRiS.rows() && "SiRiS should be square");
  //rassert(SiRiS.cols()==p && "SiRiS should have p cols");
  //rassert(isnan(SiRiS.array()).sum()==0 && "sesquare has no NaN  values");
  
}

Data_wrap::Data_wrap(const double* betahat_,const double* se_,double* sesquare_,double* q_, const size_t p_):p(p_),
												    betahat_p(betahat_),
												    se_p(se_),
												    sesquare_p(sesquare_),
												    q_p(q_)
  {

    c_mdarray betahat(betahat_p,p);
    c_mdarray se(se_p,p);
    mdarray sesquare(sesquare_p,p);
    mdarray q(q_p,p);
    
    for(size_t i=0;i<p;i++){
      sesquare(i)=se(i)*se(i); 
    }
    for(size_t i=0;i<p;i++){
      q_p[i]=betahat_p[i]/sesquare_p[i];
    }


    rassert(sesquare.cols()==1 && "sesquare should have 1 column");
    rassert(se.rows()==sesquare.rows() && "se should be same length as sesquare");
    rassert(se.size()==p && "se should have p elements");
    rassert(q.size()==p && "q should have p elements");
    rassert(betahat.size()==p && "q should have p elements");

    rassert(se.minCoeff()>0 && "se has all non-negative values");
    rassert(sesquare.minCoeff()>0 && "sesquare has all non-negative values");
    rassert(isnan(sesquare).sum()==0 && "sesquare has no NaN  values");
    rassert(isnan(q).sum()==0 && "q has no NaN  values");
    rassert(isnan(betahat).sum()==0 && "betahat has no NaN  values");


    
  }

Param_wrap::Param_wrap(double* logodds_,double* sigb_,double* sigma_beta_square_,double* s_,double* ssrat_, const Data_wrap &obj,const size_t gridsize_):
  gridsize(gridsize_),
  p(obj.p),
  logodds_p(logodds_),
  sigb_p(sigb_),
  tls(false),
  sigma_beta_square_p(sigma_beta_square_),
  s_p(s_),
  ssrat_p(ssrat_),
  logodds_t(tbbdvec(0)),
  sigma_beta_t(tbbdvec(0)),
  sigma_beta_square_t(tbbdvec(0)),
  s_t(tbbdvec(0)),
  ssrat_t(tbbdvec(0))
{

  mdarray logodds(logodds_p,gridsize);
  mdarray sigma_beta(sigb_p,gridsize);
  mdarray sigma_beta_square(sigma_beta_square_p,gridsize);
  m2darray s(s_p,p,gridsize);
  m2darray ssrat(ssrat_p,p,gridsize);
  
  Eigen::Map<Eigen::ArrayXd>  sesquare(obj.sesquare_p,p);
  
  for(size_t i=0;i<gridsize;i++){
    sigma_beta_square(i)=sigma_beta(i)*sigma_beta(i);
    s.col(i)=(sesquare*sigma_beta_square(i))/(sesquare+sigma_beta_square(i));
    ssrat.col(i)=(s.col(i)/sigma_beta_square(i)).log();
  }
  
  rassert(s.cols()==gridsize && "s should have gridsize no of columns");
  rassert(ssrat.cols()==gridsize && "ssrat should have gridsize no of columns");
  rassert(ssrat.rows()==p && "ssrat should have gridsize no of columns");
  
  rassert(sigma_beta.size()==sigma_beta_square.size() && "sigma_beta should be same length as sigma_beta_square");
  rassert(logodds.size()==gridsize && "logodds should have size of gridsize");
  
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(logodds).sum()==0 && "logodds has no NaN  values");
}


Param_wrap::Param_wrap(double* logodds_,double* sigb_,const size_t gridsize_,const size_t p_):  gridsize(gridsize_),
												p(p_),
												logodds_p(logodds_),
												sigb_p(sigb_),
												tls(true),
												sigma_beta_square_p(NULL),
												s_p(NULL),
												ssrat_p(NULL),
												logodds_t(tbbdvec(0)),
												sigma_beta_t(tbbdvec(0)),
												sigma_beta_square_t(tbbdvec(0)),
												s_t(tbbdvec(0)),
												ssrat_t(tbbdvec(0)){
  
}







Param_wrap::Param_wrap(double* sigb_,
		       double* sigma_beta_square_,
		       double* s_,
		       double* ssrat_,const Data_wrap &obj,
		       const size_t gridsize_):  gridsize(gridsize_),
						 p(obj.p),
						 logodds_p(NULL),
						 tls(false),
						 sigb_p(sigb_),
						 sigma_beta_square_p(sigma_beta_square_),
						 s_p(s_),
						 ssrat_p(ssrat_),
						 logodds_t(tbbdvec(0)),
						 sigma_beta_t(tbbdvec(0)),
						 sigma_beta_square_t(tbbdvec(0)),
						 s_t(tbbdvec(0)),
						 ssrat_t(tbbdvec(0)){
  
  mdarray sigma_beta(sigb_p,gridsize);
  mdarray sigma_beta_square(sigma_beta_square_p,gridsize);
  mdarray s(s_p,p,gridsize);
  mdarray ssrat(ssrat_p,p,gridsize);
  
  Eigen::Map<Eigen::ArrayXd>  sesquare(obj.sesquare_p,p);
  
  for(size_t i=0;i<gridsize;i++){
    sigma_beta_square(i)=sigma_beta(i)*sigma_beta(i);
    s.col(i)=(sesquare*sigma_beta_square(i))/(sesquare+sigma_beta_square(i));
    ssrat.col(i)=(s.col(i)/sigma_beta_square(i)).log();
  }
  
  rassert(s.cols()==gridsize && "s should have gridsize no of columns");
  rassert(ssrat.cols()==gridsize && "ssrat should have gridsize no of columns");
  rassert(ssrat.rows()==p && "ssrat should have gridsize no of columns");
  
  rassert(sigma_beta.size()==sigma_beta_square.size() && "sigma_beta should be same length as sigma_beta_square");
  
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  }

Param_wrap::Param_wrap(double* sigb_,const size_t gridsize_,const size_t p_):gridsize(gridsize_),
									     p(p_),
									     logodds_p(NULL),
									     sigb_p(sigb_),
									     tls(true),
									     sigma_beta_square_p(NULL),
									     s_p(NULL),
									     ssrat_p(NULL),
									     logodds_t(tbbdvec(0)),
									     sigma_beta_t(tbbdvec(0)),
									     sigma_beta_square_t(tbbdvec(0)),
									     s_t(tbbdvec(0)),
									     ssrat_t(tbbdvec(0)){

  mdarray sigma_beta(sigb_p,gridsize);
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
}


void Param_wrap::tlresize(const blocked_range<size_t> &r,const Data_wrap &obj){
  
  size_t rsize=r.end()-r.begin();
  
  logodds_t.local().resize(rsize);
  mdarray ologodds(logodds_p,gridsize);
  mdarray nlogodds(logodds_t.local().data(),rsize);
  nlogodds=ologodds.segment(r.begin(),rsize);
  //  logodds_p = logodds_t.local().data();
  
  sigma_beta_t.local().resize(rsize);
  mdarray osigma_beta(sigb_p,gridsize);
  mdarray nsigma_beta(sigma_beta_t.local().data(),rsize);
  nsigma_beta=osigma_beta.segment(r.begin(),rsize);
//  sigb_p=sigma_beta_t.local().data();
  
  sigma_beta_square_t.local().resize(rsize);
//  sigma_beta_square_p=sigma_beta_square_t.local().data();
  
  s_t.local().resize(rsize*p);
//  s_p=s_t.local().data();
  
  ssrat_t.local().resize(rsize*p);
//  ssrat_p=ssrat_t.local().data();
  
  mdarray sigma_beta(sigma_beta_t.local().data(),rsize);
  mdarray sigma_beta_square(sigma_beta_square_t.local().data(),rsize);
  m2darray s(s_t.local().data(),p,gridsize);
  m2darray ssrat(ssrat_t.local().data(),p,gridsize);
  mdarray sesquare(obj.sesquare_p,p);


  for(size_t i=r.begin();i!=r.end();i++){
    size_t ir=i-r.begin();
    //sigma_beta(ir)=sigb_p[i];
    sigma_beta_square(ir)=sigma_beta(ir)*sigma_beta(ir);
    s.col(ir)=(sesquare*sigma_beta_square(ir))/(sesquare+sigma_beta_square(ir));
    ssrat.col(ir)=(s.col(ir)/sigma_beta_square(ir)).log();
  }

  rassert((sesquare>0).all() && "sesquare is strictly positive");
  rassert((sigma_beta>0).all() && "sigma_beta is strictly positive");
  rassert((sigma_beta_square>0).all() && "sigma_beta_square is strictly positive");
  rassert(isnan(sesquare).sum()==0 && "sesquare has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(s).sum()==0 && "ssrat has no NaN  values");
  rassert(s.cols()==rsize && "s must have the right number of columns (rsize)");
  rassert(s.rows()==p && "s must have the right number of rows(p)");
  
  rassert(ssrat.cols()==rsize && "ssrat must have the right number of columns (rsize)");
  rassert(ssrat.rows()==p && "ssrat must have the right number of rows (p)");

  rassert(sigma_beta.size()==rsize && "sigma_beta must have the right number of elements (rsize)");
  rassert(sigma_beta_square.size()==rsize && "sigma_beta_square must have the right number of elements (rsize)");

  
}
  


Fit_res::Fit_res(double* lnZ_,
		 double* lnZ0_,
		 double* lnZ00_,
		 double* alpha_mean_,
		 double* pvevec_,
		 int* iternum_,
		 double* mu_mean_,
		 double* final_max_err_,
		 double* rvec_,
		 double* mtp_,
		 int* final_btnum_,const size_t p_,const size_t gridsize_):p(p_),
									   tls(false),
									   gridsize(gridsize_),lnZ_p(lnZ_),
									   lnZ0_p(lnZ0_),
									   lnZ00_p(lnZ00_),
									   alpha_mean_p(alpha_mean_),
									   pvevec_p(pvevec_),
									   iternum_p(iternum_),
									   final_max_err_p(final_max_err_),
									   mu_mean_p(mu_mean_),
									   rvec_p(rvec_),
									   final_btnum_p(final_btnum_),
									   mtp_p(mtp_),
									   lnZ_t(tbbdvec(0)),
									   lnZ0_t(tbbdvec(0)),
									   lnZ00_t(tbbdvec(0)),
									   alpha_mean_t(tbbdvec(0)),
									   mu_mean_t(tbbdvec(0)),
									   pvevec_t(tbbdvec(0)),
									   iternum_t(tbbivec(0)),
									   final_max_err_t(tbbdvec(0)),
									   final_btnum_t(tbbivec(0)),
									   rvec_t(tbbdvec(0)),
									   mtp_t(tbbdvec(0))

{


  mdarray lnZ(lnZ_p,gridsize);
  mdarray lnZ0(lnZ0_p,gridsize);
  mdarray lnZ00(lnZ00_p,gridsize);
  mdarray alpha_mean(alpha_mean_p,gridsize);
  mdarray mu_mean(mu_mean_p,gridsize);
  mdarray pvevec(pvevec_p,gridsize);
  miarray iternum(iternum_p,gridsize);

  mdarray mtp(mtp_p,gridsize);
  miarray final_btnum(final_btnum_p,gridsize);
  mdarray rvec(rvec_p,gridsize);
  mdarray final_max_err(final_max_err_p,gridsize);
  
  std::fill(iternum_p,iternum_p+gridsize,0);
  std::fill(final_btnum_p,final_btnum_p+gridsize,0);
  std::fill(alpha_mean_p,alpha_mean_p+gridsize,0);
  
  lnZ.setZero();
  lnZ0.setZero();
  lnZ00.setZero();
  mu_mean.setZero();
  rvec.setZero();
  mtp.setZero();
  
  
  rassert(iternum.sum()+final_btnum.sum()==0 && "iternum and btnum should be zeroed out");
  
  rassert(lnZ.size()==gridsize && "lnZ should be of gridsize");
  rassert(lnZ0.size()==gridsize && "lnZ0 should be of gridsize");
  rassert(lnZ00.size()==gridsize && "lnZ00 should be of gridsize");
  
  
  rassert(iternum.size()==gridsize && "iternum should be of gridsize");
  rassert(final_max_err.size()==gridsize && "final_max_err should be of gridsize");
  rassert(alpha_mean.size()==gridsize && "alpha_mean should be of gridsize");
  rassert(mu_mean.size()==gridsize && "mu_mean should be of gridsize");
  rassert(rvec.size()==gridsize && "mu_mean should be of gridsize");
  
  
  
  
}


Fit_res::Fit_res(const size_t p_,const size_t gridsize_):p(p_),gridsize(gridsize_),
							 tls(true),
							 lnZ_p(NULL),
							 lnZ0_p(NULL),
							 lnZ00_p(NULL),
							 alpha_mean_p(NULL),
							 pvevec_p(NULL),
							 iternum_p(NULL),
							 final_max_err_p(NULL),
							 mu_mean_p(NULL),
							 rvec_p(NULL),
							 final_btnum_p(NULL),
							 mtp_p(NULL),
							 lnZ_t(tbbdvec(0)),
							 lnZ0_t(tbbdvec(0)),
							 lnZ00_t(tbbdvec(0)),
							 alpha_mean_t(tbbdvec(0)),
							 mu_mean_t(tbbdvec(0)),
							 pvevec_t(tbbdvec(0)),
							 iternum_t(tbbivec(0)),
							 final_max_err_t(tbbdvec(0)),
							 final_btnum_t(tbbivec(0)),
							 rvec_t(tbbdvec(0)),
							 mtp_t(tbbdvec(0))
							 

{

  
  
  
  
}

void Fit_res::tlresize(const size_t rsize){

  
  
  lnZ_t.local().resize(rsize);
//  lnZ_p=lnZ_t.local().data();
  zero_out(lnZ_t,0);
  
  lnZ0_t.local().resize(rsize);
//  lnZ0_p=lnZ0_t.local().data();
  zero_out(lnZ0_t,0);
    
  lnZ00_t.local().resize(rsize);
//  lnZ00_p=lnZ00_t.local().data();
  zero_out(lnZ00_t,0);
  
  alpha_mean_t.local().resize(rsize);
//  alpha_mean_p=alpha_mean_t.local().data();
  zero_out(alpha_mean_t,0);
  
  mu_mean_t.local().resize(rsize);
//  mu_mean_p=mu_mean_t.local().data();
  zero_out(mu_mean_t,0);

  pvevec_t.local().resize(rsize);
//  pvevec_p=pvevec_t.local().data();
  zero_out(pvevec_t,0);
  
  iternum_t.local().resize(rsize);
//  iternum_p=iternum_t.local().data();
  zero_out(iternum_t,0);
  
  final_max_err_t.local().resize(rsize);
//  final_max_err_p=final_max_err_t.local().data();
  zero_out(final_max_err_t,0);
  
  final_btnum_t.local().resize(rsize);
//  final_btnum_p=final_btnum_t.local().data();
  zero_out(final_btnum_t,0);
  
  rvec_t.local().resize(rsize);
//  rvec_p=rvec_t.local().data();
  zero_out(rvec_t,0);
  
  mtp_t.local().resize(rsize);
//  mtp_p=mtp_t.local().data();
  zero_out(rvec_t,0);
  

  
}

rssr::rssr(const Fit_wrap &alphas_,
	   const Fit_wrap &mus_,
	   const Fit_wrap &sirisrs_,
	   const Data_wrap &datas_,
	   const Param_wrap &params_,
	   const Siris &siris_,
	   const Fit_res &results_,
	   const double tolerance_,
	   const size_t n_,
	   const size_t itermax_,
	   const std::vector<int> &forward_range_,
	   const std::vector<int> &reverse_range_): alphas(alphas_),
						    tls(alphas.tls),
						    mus(mus_),
						    sirisrs(sirisrs_),
						    results(results_),
						    datas(datas_),
						    paraams(params_),
						    siris(siris_),
						    gridsize(paraams.gridsize),
						    p(alphas.p),
						    n(n_),
						    q(datas.q_p,p),
						    logodds(paraams.logodds_p,gridsize),
						    sigma_beta(paraams.sigb_p,gridsize),
						    sigma_beta_square(paraams.sigma_beta_square_p,gridsize),
						    s(paraams.s_p,p,gridsize),
						    ssrat(paraams.ssrat_p,p,gridsize),
						    sesquare(datas.sesquare_p,p),
						    itermax(itermax_),
						    tolerance(tolerance_),
						    forward_range(forward_range_),
						    reverse_range(reverse_range_),
						    betahat(datas.betahat_p,p),
						    se(datas.se_p,p),
						    SiRiS(siris.siris_p,p,p),
						    lnZ(NULL,gridsize),
						    lnZ0(NULL,gridsize),
						    lnZ00(NULL,gridsize),
						    alpha_mean(NULL,gridsize),
						    mu_mean(NULL,gridsize),
						    pvevec(NULL,gridsize),
						    iternum(NULL,gridsize),
						    final_btnum(NULL,gridsize),
						    final_max_err(NULL,gridsize),
						    rvec(NULL,gridsize),
						    mtp(NULL,gridsize),
						    alpha(NULL,p,gridsize),
						    mu(NULL,p,gridsize),
						    SiRiSr(NULL,p,gridsize),
						    alpha0(NULL,p,gridsize),
						    mu0(NULL,p,gridsize),
						    alpha1(NULL,p,gridsize),
						    mu1(NULL,p,gridsize)

						    

{



}


void rssr::tls_init(const blocked_range<size_t> &r) const{

  size_t rstart=r.begin();
    size_t rsize= r.end()-rstart;
  if(!tls){
    rsize=gridsize;
  }
  if(tls){
    alphas.tlresize(rsize);
    mus.tlresize(rsize);
    sirisrs.tlresize(rsize);
    paraams.tlresize(r,datas);
    results.tlresize(rsize);
  }

  
  new (&alpha) m2darray(alphas.fitl.local().data(),p,rsize);
  new (&mu) m2darray(mus.fitl.local().data(),p,rsize);
  new (&alpha1) m2darray(alphas.fitl1.local().data(),p,rsize);
  new (&mu1) m2darray(mus.fitl1.local().data(),p,rsize);
  new (&alpha0) m2darray(alphas.fitl0.local().data(),p,rsize);
  new (&mu0) m2darray(mus.fitl0.local().data(),p,rsize);

  new (&SiRiSr) m2darray(sirisrs.fitl.local().data(),p,rsize);

   

  new (&lnZ) mdarray(results.lnZ_t.local().data(),rsize);
  new (&lnZ0) mdarray(results.lnZ0_t.local().data(),rsize);
  new (&lnZ00) mdarray(results.lnZ00_t.local().data(),rsize);
  new (&alpha_mean) mdarray(results.alpha_mean_t.local().data(),rsize);
  new (&mu_mean) mdarray(results.mu_mean_t.local().data(),rsize);
  new (&pvevec) mdarray(results.pvevec_t.local().data(),rsize);
  new (&iternum) miarray(results.iternum_t.local().data(),rsize);
  
  new (&mtp) mdarray(results.mtp_t.local().data(),rsize);
  new (&final_btnum) miarray(results.final_btnum_t.local().data(),rsize);
  new (&rvec) mdarray(results.rvec_t.local().data(),rsize);
  new (&final_max_err) mdarray(results.final_max_err_t.local().data(),rsize);

  new (&logodds) mdarray(paraams.logodds_t.local().data(),rsize);
  new (&sigma_beta) mdarray(paraams.sigma_beta_t.local().data(),rsize);
  new (&sigma_beta_square) mdarray(paraams.sigma_beta_square_t.local().data(),rsize);

  new (&s) m2darray(paraams.s_t.local().data(),p,rsize);
  new (&ssrat) m2darray(paraams.ssrat_t.local().data(),p,rsize);
  

  
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(logodds).sum()==0 && "logodds has no NaN  values");
    
    
  rassert(alphas.fitl.local().data()!=NULL && "alpha has to be mapped");
  rassert(mus.fitl.local().data()!=NULL && "mu has to be mapped");
  rassert(sirisrs.fitl.local().data()!=NULL && "sirisr has to be mapped");
    


  rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  rassert(isnan(alpha0).sum()==0 && "alpha0 has no NaN values");
  rassert(isnan(alpha1).sum()==0 && "alpha1 has no NaN values");
  rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
  rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");
  
 

  

  rassert(isnan(alpha).sum()==0 && "alpha has no NaN starting values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN starting values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN starting values");
}
						    

void rssr::calculate_mtp(const blocked_range <size_t > &r)const {

  for(size_t i=r.begin();i!=r.end();i++){
    mtp(i)=
      -std::sqrt(((alpha1.col(i)-alpha0.col(i)).square()).sum()+((mu1.col(i)-mu0.col(i)).square()).sum())/std::sqrt(((alpha.col(i)-2*alpha1+alpha0).square()).sum()+((mu.col(i)-2*mu1.col(i)+mu0.col(i)).square()).sum()+double_lim::epsilon());
  }
}


void rssr::calc_lnZ(const blocked_range<size_t> & r)const {
  size_t rsize=r.end()-r.begin();

  ////std::cout<<"(lnZ) rsize is :"<<rsize<<std::endl;
  rassert(lnZ.size()>=rsize && "lnZ.size() is smaller than chunk ");


  ////std::cout<<"First check that alpha is ready to go"<<std::endl;
  rassert(alpha.cols()>=rsize && "alpha has to have colno>=rsize");
  rassert(alpha.rows()==p && "alpha has to have p rows");
  rassert(alpha(0,0)!=2 && "can access alpha");
  rassert(isnan(alpha).sum()==0 && "alpha is all non NaN");

  ////std::cout<<"Now check mu"<<std::endl;
  rassert(mu.cols()>=rsize && "mu has to have colno>=rsize");
  rassert(mu.rows()==p && "mu has to have p rows");
  rassert(mu(0,0)!=2 && "can access mu");
  rassert(isnan(mu).sum()==0 && "mu is all non NaN");
  
  ////std::cout<<"Now check SiRiSr"<<std::endl;
  rassert(SiRiSr.cols()>=rsize && "SiRiSr has to have colno>=rsize");
  rassert(SiRiSr.rows()==p && "SiRiSr has to have p rows");
  rassert(SiRiSr(0,0)!=2 && "can access SiRiSr");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr is all non NaN");

  ////std::cout<<"Now check s"<<std::endl;
  rassert(s.cols()>=rsize && "s has to have colno>=rsize");
  rassert(s.rows()==p && "s has to have p rows");

  rassert(s.sum()!=0 && "can access elements of  s");
  rassert(isnan(s).sum()==0 && "s is all non NaN");
  
  ////std::cout<<"Now check ssrat"<<std::endl;
  rassert(ssrat.cols()>=rsize && "ssrat has to have colno>=rsize");
  rassert(ssrat.rows()==p && "ssrat has to have p rows");
  rassert(ssrat(0,0)!=2 && "can access ssrat");
  rassert(isnan(ssrat).sum()==0 && "ssrat is all non NaN");

  ////std::cout<<"Now check sesquare"<<std::endl;
  rassert(sesquare.size()==p && "sesquare has to have p elem");
  rassert(sesquare(0)!=2 && "can access sesquare");
  rassert(isnan(sesquare).sum()==0 && "sesquare is all non NaN");

  ////std::cout<<"Now check logodds"<<std::endl;
  rassert(logodds.size()==rsize && "logodds has to have rsize elem");
  rassert(logodds(0)!=2 && "can access logodds");
  rassert(isnan(logodds).sum()==0 && "logodds is all non NaN");

  ////std::cout<<"Now check sigma_beta"<<std::endl;
  rassert(sigma_beta.size()==rsize && "sigma_beta has to have rsize elem");
  rassert(sigma_beta(0)!=2 && "can access sigma_beta");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta is all non NaN");
  
  ////std::cout<<"about to start to calculate lnZ"<<std::endl;
    


  for(size_t j=r.begin(); j!=r.end();j++){


    

    lnZ(j)=calculate_lnZ(q,
			 alpha.col(j)*mu.col(j),
			 SiRiSr.col(j),
			 logodds(j),
			 sesquare,
			 alpha.col(j),
			 mu.col(j),
			 s.col(j),
			 sigma_beta(j));
  }
}

  

void rssr::rss_vb_iter(const blocked_range<size_t> r,bool reverse)const{
    
  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;			  
  if(reverse){

    for(const int& i: reverse_range){
      for(size_t j=r.begin();j!=r.end();j++){
	rvec(j)=alpha(i,j)*mu(i,j);
	mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	alpha(i,j)=sigmoid(logodds(j)+0.5*(ssrat(i,j)+(mu(i,j)*mu(i,j))/s(i,j)));
	SiRiSr.col(j)+=SiRiS.col(i).array()*(alpha(i,j)*mu(i,j)-rvec(j));
      }
    }

  }else{
    for(const int& i: forward_range){
      for(size_t j=r.begin();j!=r.end();j++){
	rvec(j)=alpha(i,j)*mu(i,j);
	mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	alpha(i,j)=sigmoid(logodds(j)+0.5*(ssrat(i,j)+(mu(i,j)*mu(i,j))/s(i,j)));
	SiRiSr.col(j)+=SiRiS.col(i).array()*(alpha(i,j)*mu(i,j)-rvec(j));

	
      }
	    
    }
  }
  rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
	
}

  
    
  
void rssr::squarem_step(const blocked_range<size_t> & r) const{
    
  calculate_mtp(r);
  for(size_t j=r.begin(); j!=r.end();j++){
	
    double tmtp=mtp(j);
    if(tmtp>(-1)){
      tmtp=-1;
      mtp(j) =-1;
    }else{
      if(!std::isnan(tmtp)){
	//<<"alpha_step:"<<std::endl;
	alpha.col(j)=alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*((alpha.col(j)-2*alpha1.col(j)+alpha0.col(j)));
	//<<"mu_step:"<<std::endl;
	mu.col(j)=mu0.col(j)-2*tmtp*((mu1.col(j)-mu0.col(j)))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	//<<"SiRiSr_step:"<<std::endl;
	SiRiSr.col(j)=(SiRiS.matrix()*(alpha.col(j)*mu.col(j)).matrix()).array();
      }
    }
  }
  rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
}


      
void rssr::squarem_backtrack(const blocked_range<size_t> & r)const {
      
  for(size_t j=r.begin();j!=r.end();j++){
    double tmtp=mtp(j);
    double tlnZ=lnZ(j);
    double tlnZ0=lnZ0(j);
    if(tmtp<(-1) && (tlnZ < tlnZ0)){
      size_t num_bt=0;
      while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	tmtp=0.5*(tmtp-1);
	alpha.col(j) = alpha0.col(j)-2*tmtp*(alpha1.col(j)-alpha0.col(j))+(tmtp*tmtp)*(alpha.col(j)-2*alpha1.col(j)+alpha0.col(j));
	mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	SiRiSr.col(j) = (SiRiS.matrix()*(alpha.col(j)*mu.col(j)).matrix()).array();
	rss_vb_iter(blocked_range<size_t>(j,j),false);
	rss_vb_iter(blocked_range<size_t>(j,j),true);
	      
	calc_lnZ(blocked_range<size_t>(j,j));
	num_bt=num_bt+1;
	final_btnum(j)+=1;
      }
      if(num_bt==max_backtrack){
	alpha.col(j)=alpha0.col(j);
	mu.col(j)=mu0.col(j);
	SiRiSr.col(j) = (SiRiS.matrix()*(alpha.col(j)*mu.col(j)).matrix()).array();
	calc_lnZ(blocked_range<size_t>(j,j));
	lnZ0(j)=lnZ00(j);

      }
    }
    mtp(j)=tmtp;
  }
}
    
void rssr::calc_err(const blocked_range<size_t> &r,double &max_err)const {
      
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
void rssr::compute_pve(const blocked_range<size_t> &r)const {
  // size_t rstart=r.begin();
  // size_t rsize= r.end()-rstart;
  for(size_t j=r.begin();j!=r.end();j++){
    pvevec(j)=(SiRiSr.col(j).matrix().transpose()*(alpha.col(j)*mu.col(j)).matrix());
    pvevec(j)/=(double)n;
  }
}    





  
void rssr::operator()(    const blocked_range<size_t> &qr)const {

    
  double max_err=1;
  size_t rstart=qr.begin();
  size_t rsize= qr.end()-rstart;
  blocked_range<size_t> r(0,rsize);
  tls_init(r);
  calc_lnZ(r);
  int iter=0;
  
  // Rcpp::Rcout<<"Initializing lnZ Segment "<<std::endl;    
  calc_lnZ(r);
  //Rcpp::Rcout<<"Initializing lnZ0 Segment "<<std::endl;
  lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
  while(max_err>tolerance){
    
    //Rcpp::Rcout<<"Setting lnZ00 Segment "<<std::endl;
    lnZ00.segment(rstart,rsize)=lnZ0.segment(rstart,rsize);
    //Rcpp::Rcout<<"Setting lnZ00 Segment "<<std::endl;
    lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
    //Rcpp::Rcout<<"Setting alpha0 block "<<std::endl;
    alpha0.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
    //Rcpp::Rcout<<"Setting mu0 block "<<std::endl;
    mu0.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    
    bool reverse = iter%2!=0;
    //Rcpp::Rcout<<"Performing first Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    //Rcpp::Rcout<<"Setting alpha1 block "<<std::endl;
    alpha1.block(0,rstart,p,rsize)=alpha.block(0,rstart,p,rsize);
    mu1.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    //Rcpp::Rcout<<"Performing second Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    
    rassert(isnan(alpha).sum()==0 && "alpha has no NaN  values");
    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
    
    rassert(isnan(alpha0).sum()==0 && "alpha0 has no NaN values");
    rassert(isnan(alpha1).sum()==0 && "alpha1 has no NaN values");
    


    rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
    rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");



    rassert(isnan(lnZ).sum()==0 && "lnZ has no NaN  values");
    rassert(isnan(lnZ0).sum()==0 && "lnZ0 has no NaN  values");

    //Rcpp::Rcout<<"Squarem Step "<<std::endl;
    squarem_step(r);
    //Rcpp::Rcout<<"Third Iteration "<<std::endl;
    rss_vb_iter(r,reverse);

    calc_lnZ(r);
    //Rcpp::Rcout<<"Squarem Backtrack "<<std::endl;
    squarem_backtrack(r);
    //Rcpp::Rcout<<"Error Calculation "<<std::endl;
    calc_err(r,max_err);
    if(iter>itermax){
      for(size_t  j=r.begin();j!=r.end();j++){
	final_max_err(j)=max_err;
	iternum(j)=iter;
      }
      compute_pve(r);
      break;
    }
    iter++;
  }
  //Rcpp::Rcout<<"Error Calculation "<<std::endl;
    
  for(size_t  j=r.begin();j!=r.end();j++){
    final_max_err(j)=max_err;
    iternum(j)=iter;
    alpha_mean(j)=alpha.col(j).mean();
    mu_mean(j)=mu.col(j).mean(); 
  }
  
  rassert(isnan(alpha_mean).sum()==0 && "alpha_mean is entirely finite");
  //  rassert(1==2 && "just checking that rasserts actually do something" && //std::cout<<"Something happened"<<std::endl);

  
  compute_pve(r);
}
    
  
rssr_norm::rssr_norm(const Fit_wrap &mus_,
		     const Fit_wrap &sirisrs_,
		     const Data_wrap &datas_,
		     const Param_wrap &params_,
		     const Siris &siris_,
		     const Fit_res &results_,
		     const double tolerance_, const size_t n_,
		     const size_t itermax_,
		     const std::vector<int> &forward_range_,
		     const std::vector<int> &reverse_range_
		     ): 
			mus(mus_),
			sirisrs(sirisrs_),
			results(results_),
			datas(datas_),
			paraams(params_),
			siris(siris_),
			SiRiS(siris.siris_p,p,p),
			gridsize(mus.gridsize),
			p(mus.p),
			n(n_),
			sigma_beta(NULL,gridsize),
			sigma_beta_square(NULL,gridsize),
			s(NULL,p,gridsize),
			ssrat(NULL,p,gridsize),
			itermax(itermax_),
			tolerance(tolerance_),
			forward_range(forward_range_),
			reverse_range(reverse_range_),
			tls(mus.tls),
			betahat(datas.betahat_p,p),
			se(datas.se_p,p),
			sesquare(datas.sesquare_p,p),
			q(datas.q_p,p),
			lnZ(NULL,gridsize),
			lnZ0(NULL,gridsize),
			lnZ00(NULL,gridsize),
			alpha_mean(NULL,gridsize),
			mu_mean(NULL,gridsize),
			pvevec(NULL,gridsize),
			iternum(NULL,gridsize),
			final_btnum(NULL,gridsize),
			final_max_err(NULL,gridsize),
			rvec(NULL,gridsize),
			mtp(NULL,gridsize),
			mu(NULL,p,gridsize),
			SiRiSr(NULL,p,gridsize),
			mu0(NULL,p,gridsize),
			mu1(NULL,p,gridsize)
		


{


    
       
    
}



void rssr_norm::tls_init(const blocked_range<size_t> &r) const{

  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
  if(!tls){
    rsize=gridsize;
  }
  if(tls){

    mus.tlresize(rsize);
    sirisrs.tlresize(rsize);
    paraams.tlresize(r,datas);
    results.tlresize(rsize);
  }







  new (&mu) m2darray(mus.fitl.local().data(),p,rsize);

  new (&mu1) m2darray(mus.fitl1.local().data(),p,rsize);

  new (&mu0) m2darray(mus.fitl0.local().data(),p,rsize);

  new (&SiRiSr) m2darray(sirisrs.fitl.local().data(),p,rsize);
  new (&lnZ) mdarray(results.lnZ_t.local().data(),rsize);
  new (&lnZ0) mdarray(results.lnZ0_t.local().data(),rsize);
  new (&lnZ00) mdarray(results.lnZ00_t.local().data(),rsize);
  new (&alpha_mean) mdarray(results.alpha_mean_t.local().data(),rsize);
  new (&mu_mean) mdarray(results.mu_mean_t.local().data(),rsize);
  new (&pvevec) mdarray(results.pvevec_t.local().data(),rsize);
  new (&iternum) miarray(results.iternum_t.local().data(),rsize);
  
  new (&mtp) mdarray(results.mtp_t.local().data(),rsize);
  new (&final_btnum) miarray(results.final_btnum_t.local().data(),rsize);
  new (&rvec) mdarray(results.rvec_t.local().data(),rsize);
  new (&final_max_err) mdarray(results.final_max_err_t.local().data(),rsize);


  new (&sigma_beta) mdarray(paraams.sigma_beta_t.local().data(),rsize);
  new (&sigma_beta_square) mdarray(paraams.sigma_beta_square_t.local().data(),rsize);

  new (&s) m2darray(paraams.s_t.local().data(),p,rsize);
  new (&ssrat) m2darray(paraams.ssrat_t.local().data(),p,rsize);


  

  

    
    

  rassert(mus.fitl.local().data()!=NULL && "mu has to be mapped");
  rassert(sirisrs.fitl.local().data()!=NULL && "sirisr has to be mapped");
    



  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");


  rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
  rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");
  

  rassert(isnan(mu).sum()==0 && "mu has no NaN starting values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN starting values");
}			


void rssr_norm::calculate_mtp(const blocked_range <size_t > &r)const {

  for(size_t i=r.begin();i!=r.end();i++){
    mtp(i)=-std::sqrt(((mu1.col(i)-mu0.col(i)).square()).sum())/std::sqrt(((mu.col(i)-2*mu1.col(i)+mu0.col(i)).square()).sum()+double_lim::epsilon());
  }
}


void rssr_norm::calc_lnZ(const blocked_range<size_t> & r)const {
  size_t rsize=r.end()-r.begin();

    
  rassert(lnZ.size()>=rsize && "lnZ.size() is smaller than chunk ");
  rassert(mu.cols()>=rsize && "mu has to have colno>=rsize");
  rassert(SiRiSr.cols()>=rsize && "SiRiSr has to have colno>=rsize");
  rassert(s.cols()>=rsize && "s has to have colno>=rsize");
  rassert(sigma_beta.size()>=rsize && "sigma_beta has to have size>=rsize");
  rassert(sesquare.size()==p && "sesquare has to have p elements");

  for(size_t j=r.begin(); j!=r.end();j++){
      lnZ(j)=calculate_lnZ(q,
				 mu.col(j),
				 SiRiSr.col(j),
				 sesquare,
				 mu.col(j),
				 s.col(j),
				 sigma_beta(j));
  }
}

  

void rssr_norm::rss_vb_iter(const blocked_range<size_t> r,bool reverse)const{
    
  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
    
    

  // for(size_t j=r.begin(); j!=r.end();j++){
  //        bool reverse = iternum[j]%2!=0;				  
  if(reverse){
    for(const int& i: reverse_range){
      for(size_t j=r.begin();j!=r.end();j++){
	rvec(j)=mu(i,j);
	mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	SiRiSr.col(j)+=SiRiS.col(i).array()*(mu(i,j)-rvec(j));
      }
    }

  }else{
    for(const int& i: forward_range){
      for(size_t j=r.begin();j!=r.end();j++){
	rvec(j)=mu(i,j);
	mu(i,j)=s(i,j)*(q(i)+rvec(j)/sesquare(i)-SiRiSr(i,j));
	SiRiSr.col(j)+=SiRiS.col(i).array()*(mu(i,j)-rvec(j));
      }	    
    }
  }

  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
	
}

    
    
  
void rssr_norm::squarem_step(const blocked_range<size_t> & r) const{
    
  calculate_mtp(r);
  for(size_t j=r.begin(); j!=r.end();j++){
	
    double tmtp=mtp(j);
    if(tmtp>(-1)){
      tmtp=-1;
      mtp(j) =-1;
    }else{
      if(!std::isnan(tmtp)){
	mu.col(j)=mu0.col(j)-2*tmtp*((mu1.col(j)-mu0.col(j)))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	//Rcpp::Rcout<<"SiRiSr_step:"<<std::endl;
	SiRiSr.col(j)=(SiRiS.matrix()*(mu.col(j)).matrix()).array();
      }
    }
  }

  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
}


      
void rssr_norm::squarem_backtrack(const blocked_range<size_t> & r)const {
      
  for(size_t j=r.begin();j!=r.end();j++){
    double tmtp=mtp(j);
    double tlnZ=lnZ(j);
    double tlnZ0=lnZ0(j);
    if(tmtp<(-1) && (tlnZ < tlnZ0)){
      size_t num_bt=0;
      while((tlnZ<tlnZ0) && (num_bt < max_backtrack)){
	tmtp=0.5*(tmtp-1);
	mu.col(j) = mu0.col(j)-2*tmtp*(mu1.col(j)-mu0.col(j))+(tmtp*tmtp)*(mu.col(j)-2*mu1.col(j)+mu0.col(j));
	SiRiSr.col(j) = (SiRiS.matrix()*(mu.col(j)).matrix()).array();
	rss_vb_iter(blocked_range<size_t>(j,j),false);
	rss_vb_iter(blocked_range<size_t>(j,j),true);
	calc_lnZ(blocked_range<size_t>(j,j));
	num_bt=num_bt+1;
	final_btnum(j)++;
      }
      if(num_bt==max_backtrack){
	mu.col(j)=mu0.col(j);
	SiRiSr.col(j) = (SiRiS.matrix()*(mu.col(j)).matrix()).array();
	calc_lnZ(blocked_range<size_t>(j,j));
	lnZ0(j)=lnZ00(j);
      }
    }
    mtp(j)=tmtp;
  }
}
    
void rssr_norm::calc_err(const blocked_range<size_t> &r,double &max_err)const {
  Eigen::ArrayXd as(p);
  as.setOnes();
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
      double tmax_err=find_maxerr(as,as,mu.col(j),mu0.col(j));
      if(tmax_err>max_err){
	max_err=tmax_err;
      }
    }
  }
}

void rssr_norm::compute_pve(const blocked_range<size_t> &r)const {

  for(size_t j=r.begin();j!=r.end();j++){
    pvevec(j)=(SiRiSr.col(j).matrix().transpose()*(mu.col(j)).matrix());
    pvevec(j)/=(double)n;
  }
} 



void rssr_norm::operator()(    const blocked_range<size_t> &qr)const {
    

  double max_err=1;
  size_t rstart=qr.begin();
  size_t rsize= qr.end()-rstart;
   

  tls_init(qr);
  blocked_range<size_t> r(0,rsize);
  calc_lnZ(r);
  int iter=0;

  
  ////Rcpp::Rcout<<"Initializing lnZ Segment "<<std::endl;    
  //  calc_lnZ(r);
  //<<"Initializing lnZ0 Segment "<<std::endl;
  lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
  while(max_err>tolerance){

    //<<"Setting lnZ00 Segment "<<std::endl;
    lnZ00.segment(rstart,rsize)=lnZ0.segment(rstart,rsize);
    //<<"Setting lnZ00 Segment "<<std::endl;
    lnZ0.segment(rstart,rsize)=lnZ.segment(rstart,rsize);
    //<<"Setting alpha0 block "<<std::endl;
    //<<"Setting mu0 block "<<std::endl;
    mu0.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);

    bool reverse = iter%2!=0;
    //<<"Performing first Iteration "<<std::endl;
    rss_vb_iter(r,reverse);
    //<<"Setting alpha1 block "<<std::endl;

    mu1.block(0,rstart,p,rsize)=mu.block(0,rstart,p,rsize);
    //<<"Performing second Iteration "<<std::endl;
    rss_vb_iter(r,reverse);


    rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
    rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");



    rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
    rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");



    rassert(isnan(lnZ).sum()==0 && "lnZ has no NaN  values");
    rassert(isnan(lnZ0).sum()==0 && "lnZ0 has no NaN  values");


    squarem_step(r);
    //<<"Third Iteration "<<std::endl;
    rss_vb_iter(r,reverse);

    calc_lnZ(r);
    //<<"Squarem Backtrack "<<std::endl;
    squarem_backtrack(r);
    //<<"Error Calculation "<<std::endl;
    calc_err(r,max_err);
    if(iter>itermax){
      for(size_t  j=r.begin();j!=r.end();j++){
	final_max_err(j)=max_err;
	iternum(j)=iter;
      }
      compute_pve(r);
      break;
    }
    iter++;
  }
    
  for(size_t  j=r.begin();j!=r.end();j++){
    final_max_err(j)=max_err;
    iternum(j)=iter;
    alpha_mean(j)=1;
    mu_mean(j)=mu.col(j).mean();
  }
  rassert(isnan(alpha_mean).sum()==0 && "alpha_mean is entirely finite");
  compute_pve(r);
}
