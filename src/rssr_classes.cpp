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


Fit_wrap::Fit_wrap(const double* fit_,const double* fit0_,const double* fit1_,const size_t p_, const size_t gridsize_): p(p_),
												      gridsize(gridsize_),
												      keepold(true),
												      ofit_p(fit_),
												      fit_t(tbbdvec(0)),
												      fit0_t(tbbdvec(0)),
												      fit1_t(tbbdvec(0)),
															fit_m(m2darray(fit_t.local().data(),p,gridsize)),
															fit0_m(m2darray(fit0_t.local().data(),p,gridsize)),
															fit1_m(m2darray(fit1_t.local().data(),p,gridsize)){
#ifndef NDEBUG
  mdarray fit(ofit_p,p);
  rassert(fit.size()==p &&"fit should have size of number of SNPs (p)");
  rassert(isnan(fit).sum()==0 && "fit has no NaN  values");
#endif
    
  }

Fit_wrap::Fit_wrap(const double* fit_,const size_t p_, const size_t gridsize_): p(p_),
									  gridsize(gridsize_),
									  keepold(false),
									  ofit_p(fit_),
									  fit_t(tbbdvec(0)),
									  fit0_t(tbbdvec(0)),
									  fit1_t(tbbdvec(0)),
										fit_m(m2darray(fit_t.local().data(),p,gridsize)),
										fit0_m(m2darray(NULL,p,gridsize)),
										fit1_m(m2darray(NULL,p,gridsize)){
#ifndef NDEBUG
  mdarray fit(ofit_p,p);
  rassert(fit.size()==p &&"fit should have size of number of SNPs (p)");
  rassert(isnan(fit).sum()==0 && "fit has no NaN  values");
#endif
    
}


void Fit_wrap::tlresize(const size_t rsize){
  

  fit_t.local().reserve(rsize*p);
  c_mdarray ofit(ofit_p,p);
  //  m2darray nfit(fit_t.local().data(),p,rsize);
  new (&fit_m.local()) m2darray(fit_t.local().data(),p,rsize);
  for(int i=0; i<rsize;i++){
    fit_m.local().col(i)=ofit;
  }
  rassert(fit_m.local().cols()==rsize &&"fit should have same num of cols as rsize");
  rassert(isnan(fit_m.local()).sum()==0 && "fit has no NaN  values");

  //  fit_p=fitl.local().data();
  
  if(keepold){
    //    Rcpp::Rcout<<"resizing fit0"<<std::endl;
    fit0_t.local().resize(rsize*p);
    new (&fit0_m.local()) m2darray(fit0_t.local().data(),p,rsize);
    
    fit1_t.local().resize(rsize*p);
    new (&fit1_m.local()) m2darray(fit1_t.local().data(),p,rsize);
    
    
  }
}




Data_wrap::Data_wrap(const double* betahat_,const double* se_,const double* sesquare_,const double* q_,const double *siris_, const size_t p_):p(p_),
																	      betahat_p(betahat_),
																	      se_p(se_),
																	      siris_p(siris_),
																	      sesquare_p(sesquare_),
																	      q_p(q_){
  c_mdarray betahat (betahat_p,p);
  c_mdarray se (se_p,p);
  c_m2darray siris (siris_p,p,p);
  c_mdarray sesquare (sesquare_p,p);
  c_mdarray q(q_p,p);

		      
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

Param_wrap::Param_wrap(const double* logodds_,const double* sigb_,const double* sigma_beta_square_,const size_t gridsize_,const size_t p_):gridsize(gridsize_),
															p(p_),
															ologodds_p(logodds_),
															osigma_beta_p(sigb_),
															osigma_beta_square_p(sigma_beta_square_),
															logodds_t(tbbdvec(0)),
															sigma_beta_t(tbbdvec(0)),
															sigma_beta_square_t(tbbdvec(0)),
															s_t(tbbdvec(0)),
															ssrat_t(tbbdvec(0)),
																	   logodds_m(mdarray(logodds_t.local().data(),gridsize)),
																	   sigma_beta_m(mdarray(sigma_beta_t.local().data(),gridsize)),
																	   sigma_beta_square_m(mdarray(sigma_beta_square_t.local().data(),gridsize)),
																	   s_m(m2darray(s_t.local().data(),p,gridsize)),
																	   ssrat_m(m2darray(ssrat_t.local().data(),p,gridsize)){
  

  c_mdarray logodds(ologodds_p,gridsize);
  c_mdarray sigma_beta(osigma_beta_p,gridsize);
  c_mdarray sigma_beta_square(osigma_beta_square_p,gridsize);
  
  
  
  
  rassert(sigma_beta.size()==sigma_beta_square.size() && "sigma_beta should be same length as sigma_beta_square");
  rassert(logodds.size()==gridsize && "logodds should have size of gridsize");
  
  
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  // rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(logodds).sum()==0 && "logodds has no NaN  values");
}


Param_wrap::Param_wrap(const double* sigb_,const double* sigma_beta_square_,const size_t gridsize_,const size_t p_):gridsize(gridsize_),
														    p(p_),
														    ologodds_p(NULL),
														    osigma_beta_p(sigb_),
														    osigma_beta_square_p(sigma_beta_square_),
														    logodds_t(tbbdvec(0)),
														    sigma_beta_t(tbbdvec(0)),
														    sigma_beta_square_t(tbbdvec(0)),
														    s_t(tbbdvec(0)),
														    ssrat_t(tbbdvec(0)),
														    logodds_m(mdarray(logodds_t.local().data(),gridsize)),
														    sigma_beta_m(mdarray(sigma_beta_t.local().data(),gridsize)),
														    sigma_beta_square_m(mdarray(sigma_beta_square_t.local().data(),gridsize)),
														    s_m(m2darray(s_t.local().data(),p,gridsize)),
														    ssrat_m(m2darray(ssrat_t.local().data(),p,gridsize)){
  
  
  
  c_mdarray sigma_beta(osigma_beta_p,gridsize);
  c_mdarray sigma_beta_square(osigma_beta_square_p,gridsize);
  

  
  rassert(sigma_beta.size()==sigma_beta_square.size() && "sigma_beta should be same length as sigma_beta_square");

  

  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  //rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
 }





void Param_wrap::tlresize(const blocked_range<size_t> &r,const Data_wrap &obj){
  
  size_t rsize=r.end()-r.begin();
  
  logodds_t.local().resize(rsize);
  if(ologodds_p!=NULL){
    c_mdarray ologodds(ologodds_p,gridsize);

    new (&logodds_m.local()) mdarray(logodds_t.local().data(),rsize);
    //        mdarray nlogodds(logodds_t.local().data(),rsize);
    logodds_m.local()=ologodds.segment(r.begin(),rsize);
  }
  //  logodds_p = logodds_t.local().data();
  
  sigma_beta_t.local().resize(rsize);
  c_mdarray osigma_beta(osigma_beta_p,gridsize);


  new (&sigma_beta_m.local()) mdarray(sigma_beta_t.local().data(),rsize);
  sigma_beta_m.local()=osigma_beta.segment(r.begin(),rsize);
 std::cout<<"rstart:"<<r.begin()<<" rsize:"<<rsize<<endl;
 std::cout<<"Sigma_beta_m.local():"<<sigma_beta_m.local()<<std::endl;
 std::cout<<"osigma_beta.segment:"<<osigma_beta.segment(r.begin(),rsize)<<std::endl;
  
  //  mdarray nsigma_beta(sigma_beta_t.local().data(),rsize);
  //  sigb_p=sigma_beta_t.local().data();

  c_mdarray osigma_beta_square(osigma_beta_square_p,gridsize);
  sigma_beta_square_t.local().resize(rsize);
  new (&sigma_beta_square_m.local()) mdarray(sigma_beta_square_t.local().data(),rsize);
  sigma_beta_square_m.local()=osigma_beta_square.segment(r.begin(),rsize);
  //  sigma_beta_square_p=sigma_beta_square_t.local().data();
  
  s_t.local().resize(rsize*p);
  new (&s_m.local()) m2darray(s_t.local().data(),p,rsize);

  //  s_p=s_t.local().data();
  
  ssrat_t.local().resize(rsize*p);
  new (&ssrat_m.local()) m2darray(ssrat_t.local().data(),p,rsize);

  
  //  mdarray sigma_beta(sigma_beta_t.local().data(),rsize);
  //  mdarray sigma_beta_square(sigma_beta_square_t.local().data(),rsize);
  //  m2darray s(s_t.local().data(),p,gridsize);
  //  m2darray ssrat(ssrat_t.local().data(),p,gridsize);
  c_mdarray sesquare(obj.sesquare_p,p);


  for(size_t i=r.begin();i!=r.end();i++){
    size_t ir=i-r.begin();
    //sigma_beta(ir)=sigb_p[i];
    //    sigma_beta_square_m.local().(ir)=sigma_beta_m.local().(ir)*sigma_beta_m.local().(ir);
    s_m.local().col(ir)=(sesquare*sigma_beta_square_m.local().coeff(ir))/(sesquare+sigma_beta_square_m.local().coeff(ir));
    ssrat_m.local().col(ir)=(s_m.local().col(ir)/sigma_beta_square_m.local().coeff(ir)).log();
  }

  rassert((sesquare>0).all() && "sesquare is strictly positive");
  rassert((sigma_beta_m.local()>0).all() && "sigma_beta is strictly positive");
  rassert((sigma_beta_square_m.local()>0).all() && "sigma_beta_square is strictly positive");
  rassert(isnan(sesquare).sum()==0 && "sesquare has no NaN  values");
  rassert(isnan(sigma_beta_square_m.local()).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(sigma_beta_m.local()).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(ssrat_m.local()).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(s_m.local()).sum()==0 && "ssrat has no NaN  values");
  rassert(s_m.local().cols()==rsize && "s must have the right number of columns (rsize)");
  rassert(s_m.local().rows()==p && "s must have the right number of rows(p)");
  
  rassert(ssrat_m.local().cols()==rsize && "ssrat must have the right number of columns (rsize)");
  rassert(ssrat_m.local().rows()==p && "ssrat must have the right number of rows (p)");

  rassert(sigma_beta_m.local().size()==rsize && "sigma_beta must have the right number of elements (rsize)");
  rassert(sigma_beta_square_m.local().size()==rsize && "sigma_beta_square must have the right number of elements (rsize)");
}

  
Fit_res::Fit_res(const size_t p_,
		 const size_t gridsize_):p(p_),
					 gridsize(gridsize_),
					 lnZ_p(NULL),
					 alpha_mean_p(NULL),
					 pvevec_p(NULL),
					 iternum_p(NULL),
					 final_max_err_p(NULL),
					 mu_mean_p(NULL),
					 final_btnum_p(NULL),
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
					 mtp_t(tbbdvec(0)),
					 lnZ_m(mdarray(NULL,gridsize)),
					 lnZ0_m(mdarray(NULL,gridsize)),
					 lnZ00_m(mdarray(NULL,gridsize)),
					 alpha_mean_m(mdarray(NULL,gridsize)),
					 mu_mean_m(mdarray(NULL,gridsize)),
					 pvevec_m(mdarray(NULL,gridsize)),
					 iternum_m(miarray(NULL,gridsize)),
					 final_max_err_m(mdarray(NULL,gridsize)),
					 final_btnum_m(miarray(NULL,gridsize)),
					 rvec_m(mdarray(NULL,gridsize)),
					 mtp_m(mdarray(NULL,gridsize))
{
  


  
  
  
}

Fit_res::Fit_res(double* lnZ_,
	double* alpha_mean_,
	double* pvevec_,
	int* iternum_,
	double* mu_mean_,
	double* final_max_err_,
	int* final_btnum_,const size_t p_,const size_t gridsize_):p(p_),
								  gridsize(gridsize_),
								  lnZ_p(lnZ_),
								  alpha_mean_p(alpha_mean_),
								  pvevec_p(pvevec_),
								  iternum_p(iternum_),
								  final_max_err_p(final_max_err_),
								  mu_mean_p(mu_mean_),
								  final_btnum_p(final_btnum_),
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
								  mtp_t(tbbdvec(0)),
								  lnZ_m(mdarray(NULL,gridsize)),
								  lnZ0_m(mdarray(NULL,gridsize)),
								  lnZ00_m(mdarray(NULL,gridsize)),
								  alpha_mean_m(mdarray(NULL,gridsize)),
								  mu_mean_m(mdarray(NULL,gridsize)),
								  pvevec_m(mdarray(NULL,gridsize)),
								  iternum_m(miarray(NULL,gridsize)),
								  final_max_err_m(mdarray(NULL,gridsize)),
								  final_btnum_m(miarray(NULL,gridsize)),
								  rvec_m(mdarray(NULL,gridsize)),
								  mtp_m(mdarray(NULL,gridsize))
{
  


  
  
  
}



  
  
  
  


void Fit_res::tlresize(const size_t rsize){

  
  
  lnZ_t.local().resize(rsize);
  //  lnZ_p=lnZ_t.local().data();
  new (&lnZ_m.local()) mdarray(lnZ_t.local().data(),rsize);
  
  lnZ0_t.local().resize(rsize);
//  lnZ0_p=lnZ0_t.local().data();
  new (&lnZ0_m.local()) mdarray(lnZ0_t.local().data(),rsize);
    
  lnZ00_t.local().resize(rsize);
//  lnZ00_p=lnZ00_t.local().data();
  new (&lnZ00_m.local()) mdarray(lnZ00_t.local().data(),rsize);
  
  alpha_mean_t.local().resize(rsize);
//  alpha_mean_p=alpha_mean_t.local().data();
  new (&alpha_mean_m.local()) mdarray(alpha_mean_t.local().data(),rsize);
  
  mu_mean_t.local().resize(rsize);
//  mu_mean_p=mu_mean_t.local().data();
  new (&mu_mean_m.local()) mdarray(mu_mean_t.local().data(),rsize);

  pvevec_t.local().resize(rsize);
//  pvevec_p=pvevec_t.local().data();
  new (&pvevec_m.local()) mdarray(pvevec_t.local().data(),rsize);
  
  iternum_t.local().resize(rsize);
//  iternum_p=iternum_t.local().data();
  new (&iternum_m.local()) miarray(iternum_t.local().data(),rsize);
  
  final_max_err_t.local().resize(rsize);
//  final_max_err_p=final_max_err_t.local().data();
  new (&final_max_err_m.local()) mdarray(final_max_err_t.local().data(),rsize);
  
  final_btnum_t.local().resize(rsize);
//  final_btnum_p=final_btnum_t.local().data();
  new (&final_btnum_m.local()) miarray(final_btnum_t.local().data(),rsize);
  
  rvec_t.local().resize(rsize);
//  rvec_p=rvec_t.local().data();
  new (&rvec_m.local()) mdarray(rvec_t.local().data(),rsize);
  
  mtp_t.local().resize(rsize);
  //  mtp_p=mtp_t.local().data();
  new (&mtp_m.local()) mdarray(mtp_t.local().data(),rsize);
 
  
}

void Fit_res::write_res(const blocked_range<size_t> & r){

  if(lnZ_p==NULL){
    Rcpp::stop("lnZ_p not allocated!");
  }
  
  for(size_t i=r.begin();i!=r.end();i++){
    size_t ir=i-r.begin();
    lnZ_p[i]=lnZ_m.local().coeff(ir);
    alpha_mean_p[i]=alpha_mean_m.local().coeff(ir);
    mu_mean_p[i]=mu_mean_m.local().coeff(ir);
    pvevec_p[i]=pvevec_m.local().coeff(ir);
    final_max_err_p[i]=final_max_err_m.local().coeff(ir);
    iternum_p[i]=iternum_m.local().coeff(ir);
    final_btnum_p[i]=final_btnum_m.local().coeff(ir);
  }
}
    
    
    
  
  

rssr::rssr(Fit_wrap *alphas_,
	   Fit_wrap *mus_,
	   Fit_wrap *sirisrs_,
	   const Data_wrap *datas_,
	   Param_wrap *params_,
	   Fit_res *results_,
	   const double tolerance_,
	   const size_t n_,
	   const size_t itermax_,
	   const std::vector<int> &forward_range_,
	   const std::vector<int> &reverse_range_): alphas(alphas_),
						    tls(true),
						    mus(mus_),
						    sirisrs(sirisrs_),
						    results(results_),
						    datas(datas_),
						    paraams(params_),
						    gridsize(paraams->gridsize),
						    p(alphas->p),
						    n(n_),
						    q(datas->q_p,p),
						    logodds(paraams->logodds_m.local()),
						    sigma_beta(paraams->sigma_beta_m.local()),
						    sigma_beta_square(paraams->sigma_beta_square_m.local()),
						    s(paraams->s_m.local()),
						    ssrat(paraams->ssrat_m.local()),
						    sesquare(datas->sesquare_p,p),
						    itermax(itermax_),
						    tolerance(tolerance_),
						    forward_range(forward_range_),
						    reverse_range(reverse_range_),
						    betahat(datas->betahat_p,p),
						    se(datas->se_p,p),
						    SiRiS(datas->siris_p,p,p),
						    lnZ(results->lnZ_m.local()),
						    lnZ0(results->lnZ0_m.local()),
						    lnZ00(results->lnZ00_m.local()),
						    alpha_mean(results->alpha_mean_m.local()),
						    mu_mean(results->mu_mean_m.local()),
						    pvevec(results->pvevec_m.local()),
						    iternum(results->iternum_m.local()),
						    final_btnum(results->final_btnum_m.local()),
						    final_max_err(results->final_max_err_m.local()),
						    rvec(results->rvec_m.local()),
						    mtp(results->mtp_m.local()),
						    alpha(alphas->fit_m.local()),
						    mu(mus->fit_m.local()),
						    SiRiSr(sirisrs->fit_m.local()),
						    alpha0(alphas->fit0_m.local()),
						    mu0(mus->fit0_m.local()),
						    alpha1(alphas->fit1_m.local()),
						    mu1(mus->fit1_m.local())

						    

{



}


void rssr::tls_init(const blocked_range<size_t> &r) const{

  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;
  if(!tls){
    rsize=gridsize;
  }
  if(tls){
    alphas->tlresize(rsize);
    mus->tlresize(rsize);
    sirisrs->tlresize(rsize);
    paraams->tlresize(r,*datas);
    results->tlresize(rsize);
  }

  
  alpha = alphas->fit_m.local();
  mu = mus->fit_m.local();
  alpha1 = alphas->fit1_m.local();
  mu1 = mus->fit1_m.local();
  alpha0 = alphas->fit0_m.local();
  mu0 = mus->fit0_m.local();

  SiRiSr = sirisrs->fit_m.local();

   

  lnZ = results->lnZ_m.local();
  lnZ0 = results->lnZ0_m.local();
  lnZ00 = results->lnZ00_m.local();
  alpha_mean = results->alpha_mean_m.local();
  mu_mean = results->mu_mean_m.local();
  pvevec = results->pvevec_m.local();
  iternum = results->iternum_m.local();
  
  mtp = results->mtp_m.local();
  final_btnum = results->final_btnum_m.local();
  rvec = results->rvec_m.local();
  final_max_err = results->final_max_err_m.local();

  logodds = paraams->logodds_m.local();
  sigma_beta = paraams->sigma_beta_m.local();
  sigma_beta_square = paraams->sigma_beta_square_m.local();

  s = paraams->s_m.local();
  ssrat = paraams->ssrat_m.local();
  

  
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");
  rassert(isnan(logodds).sum()==0 && "logodds has no NaN  values");
    

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
      results->write_res(qr);
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
  results->write_res(qr);
}
    
  
rssr_norm::rssr_norm(Fit_wrap *mus_,
		     Fit_wrap *sirisrs_,
		     const Data_wrap *datas_,
		     Param_wrap *params_,
		     Fit_res *results_,
		     const double tolerance_, const size_t n_,
		     const size_t itermax_,
		     const std::vector<int> &forward_range_,
		     const std::vector<int> &reverse_range_
		     ): tls(true),
			mus(mus_),
			sirisrs(sirisrs_),
			results(results_),
			datas(datas_),
			paraams(params_),
			gridsize(paraams->gridsize),
			p(mus->p),
			n(n_),
			q(datas->q_p,p),
			sigma_beta(paraams->sigma_beta_m.local()),
			sigma_beta_square(paraams->sigma_beta_square_m.local()),
			s(paraams->s_m.local()),
			ssrat(paraams->ssrat_m.local()),
			sesquare(datas->sesquare_p,p),
			itermax(itermax_),
			tolerance(tolerance_),
			forward_range(forward_range_),
			reverse_range(reverse_range_),
			betahat(datas->betahat_p,p),
			se(datas->se_p,p),
			SiRiS(datas->siris_p,p,p),
			lnZ(results->lnZ_m.local()),
			lnZ0(results->lnZ0_m.local()),
			lnZ00(results->lnZ00_m.local()),
			alpha_mean(results->alpha_mean_m.local()),
			mu_mean(results->mu_mean_m.local()),
			pvevec(results->pvevec_m.local()),
			iternum(results->iternum_m.local()),
			final_btnum(results->final_btnum_m.local()),
			final_max_err(results->final_max_err_m.local()),
			rvec(results->rvec_m.local()),
			mtp(results->mtp_m.local()),
			mu(mus->fit_m.local()),
			SiRiSr(sirisrs->fit_m.local()),
			mu0(mus->fit0_m.local()),
			mu1(mus->fit1_m.local()){


    
       
    
}



void rssr_norm::tls_init(const blocked_range<size_t> &r) const{



  size_t rstart=r.begin();
  size_t rsize= r.end()-rstart;

  //    rsize=gridsize;
  mus->tlresize(rsize);
  sirisrs->tlresize(rsize);
  paraams->tlresize(r,*datas);
  results->tlresize(rsize);

  
  
  mu = mus->fit_m.local();

  mu1 = mus->fit1_m.local();
 
  mu0 = mus->fit0_m.local();

  SiRiSr = sirisrs->fit_m.local();

   

  lnZ = results->lnZ_m.local();
  lnZ0 = results->lnZ0_m.local();
  lnZ00 = results->lnZ00_m.local();
  alpha_mean = results->alpha_mean_m.local();
  mu_mean = results->mu_mean_m.local();
  pvevec = results->pvevec_m.local();
  iternum = results->iternum_m.local();
  
  mtp = results->mtp_m.local();
  final_btnum = results->final_btnum_m.local();
  rvec = results->rvec_m.local();
  final_max_err = results->final_max_err_m.local();


  sigma_beta = paraams->sigma_beta_m.local();
  sigma_beta_square = paraams->sigma_beta_square_m.local();

  s = paraams->s_m.local();
  ssrat = paraams->ssrat_m.local();
  

  
  rassert(isnan(s).sum()==0 && "s has no NaN  values");
  rassert(isnan(sigma_beta).sum()==0 && "sigma_beta has no NaN  values");
  rassert(isnan(sigma_beta_square).sum()==0 && "sigma_beta_square has no NaN  values");
  rassert(isnan(ssrat).sum()==0 && "ssrat has no NaN  values");

    

  rassert(isnan(mu).sum()==0 && "mu has no NaN  values");
  rassert(isnan(SiRiSr).sum()==0 && "SiRiSr has no NaN  values");
  rassert(isnan(mu0).sum()==0 && "mu0 has no NaN  values");
  rassert(isnan(mu1).sum()==0 && "mu1 has no NaN  values");
 






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
