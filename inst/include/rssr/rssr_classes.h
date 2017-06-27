#ifndef INCLUDE_TBB
#define INCLUDE_TBB
#include <RcppEigen.h>
#include "rssr_types.h"


using namespace Eigen;

#include <cmath>
#include <RcppParallel.h>
#include <algorithm>

using namespace std;



#include <stdio.h>



using namespace RcppParallel;
using namespace tbb;






class Fit_wrap{
 public:  
  const size_t p;
  const size_t gridsize;
  const bool keepold;
  const double const *ofit_p;

  fitdtype fit_t;
  fitdtype fit1_t;
  fitdtype fit0_t;

  tl2darray fit_m;
  tl2darray fit0_m;
  tl2darray fit1_m;
  
  Fit_wrap(const double* fit_,const double* fit0_,const double* fit1_,const size_t p_, const size_t gridsize_);
  Fit_wrap(const double* fit_,const size_t p_,const size_t gridsize);
  void tlresize(const size_t rsize);
};

class Data_wrap{
 public:
  const size_t p;
  const double *const betahat_p;
  const double *const se_p;
  const double *const siris_p;
  const double* const sesquare_p;
  const double* const q_p;
  Data_wrap(const double* betahat_,const double* se_,const double* sesquare_,const double* q_,const double* siris_, const size_t p_);
												   
};




class Param_wrap{
 public:
  const size_t gridsize;
  const size_t p;
  const double *const ologodds_p;
  const double *const osigma_beta_p;
  const double *const osigma_beta_square_p;

  fitdtype logodds_t;
  fitdtype sigma_beta_t;
  fitdtype sigma_beta_square_t;
  fitdtype s_t;
  fitdtype ssrat_t;

  tldarray logodds_m;
  tldarray sigma_beta_m;
  tldarray sigma_beta_square_m;
  tl2darray s_m;
  tl2darray ssrat_m;

  
  Param_wrap(const double* logodds_,const double* sigb_,const double *sigma_beta_square_, const size_t gridsize_,const size_t p_);
  Param_wrap(const double* sigb_,const double* sigma_beta_square,const size_t gridsize_,const size_t p_);
  void tlresize(const blocked_range<size_t> &r,const Data_wrap &obj);
};


class Fit_res{
 public:
  const size_t gridsize;
  const size_t p;

  double *lnZ_p;
  double *alpha_mean_p;
  double *mu_mean_p;
  double *pvevec_p;
  double *final_max_err_p;
  
  int    *iternum_p;
  int    *final_btnum_p;



  fitdtype lnZ_t;
  fitdtype lnZ0_t;
  fitdtype lnZ00_t;
  fitdtype alpha_mean_t;
  fitdtype mu_mean_t;
  fitdtype pvevec_t;
  fititype iternum_t;
  fitdtype final_max_err_t;
  fititype final_btnum_t;
  fitdtype rvec_t;
  fitdtype mtp_t;

  
  tldarray lnZ_m;
  tldarray lnZ0_m;
  tldarray lnZ00_m;
  tldarray alpha_mean_m;
  tldarray mu_mean_m;
  tldarray pvevec_m;
  tliarray iternum_m;
  tldarray final_max_err_m;
  tliarray final_btnum_m;
  tldarray rvec_m;
  tldarray mtp_m;


  

  
  Fit_res(double* lnZ_,
	  double* alpha_mean_,
	  double* pvevec_,
	  int* iternum_,
	  double* mu_mean_,
	  double* final_max_err_,
	  int* final_btnum_,const size_t p_,const size_t gridsize_);
  Fit_res(const size_t p_,const size_t gridsize_);
  void write_res(const blocked_range<size_t> & r);
  void tlresize (const size_t rsize);


};



class rssr_norm{
 public:
  
  mutable Fit_wrap *mus;
  mutable Fit_wrap *sirisrs;
  mutable Fit_res *results;
  
  const bool use_lnztol=true;
  const bool tls;
  const size_t p;
  const size_t n;


  const int itermax;
  const double tolerance;

  mutable Param_wrap *paraams;
  const Data_wrap *const datas;
  const size_t gridsize;


  
  tl2darray::reference mu;
  tl2darray::reference SiRiSr;

  
  tl2darray::reference mu0;

  
  tl2darray::reference mu1;

  const c_mdarray betahat;
  const c_mdarray se;
  const c_mdarray sesquare;
  const c_mdarray q;
  const c_m2darray SiRiS;

  
  
  tldarray::reference sigma_beta;
  tldarray::reference sigma_beta_square;
  tl2darray::reference s;
  tl2darray::reference ssrat;
  
  tldarray::reference lnZ;
  tldarray::reference lnZ0;
  tldarray::reference lnZ00;
  tldarray::reference alpha_mean;
  tldarray::reference mu_mean;
  tldarray::reference pvevec;
  
  tliarray::reference iternum;
  tliarray::reference final_btnum;
  
  tldarray::reference final_max_err;
  tldarray::reference rvec;
  tldarray::reference mtp;
  
  
  std::vector<int> forward_range;
  std::vector<int> reverse_range;

  const size_t max_backtrack=12;
  

  rssr_norm(Fit_wrap *mus_,
	    Fit_wrap *sirisrs_,
	    const Data_wrap *datas_,
	    Param_wrap *params_,
	    Fit_res *results_,
	    const double tolerance_, const size_t n_,
	    const size_t itermax_,
	    const std::vector<int> &forward_range_,
	    const std::vector<int> &reverse_range_
	    );
  void calc_lnZ(const blocked_range<size_t> & r)const;
  void calculate_mtp(const blocked_range <size_t > &r)const;
  void rss_vb_iter(const blocked_range<size_t> r,bool reverse)const;
  void squarem_step(const blocked_range<size_t> & r)const;
  void squarem_backtrack(const blocked_range<size_t> & r)const;
  void calc_err(const blocked_range<size_t> &r,double &max_err)const;
  void compute_pve(const blocked_range<size_t> &r)const;
  void operator()(    const blocked_range<size_t> &r)const;
  void tls_init(const blocked_range<size_t> &r)const;
  
  
};



class rssr{
 public:
  
  mutable Fit_wrap *alphas;
  mutable Fit_wrap *mus;
  mutable Fit_wrap *sirisrs;
  mutable Fit_res *results;
  
  const bool use_lnztol=true;
  const bool tls;
  const size_t p;
  const size_t n;


  const int itermax;
  const double tolerance;

  mutable Param_wrap *paraams;
  const Data_wrap *const datas;
  size_t gridsize;

  tl2darray::reference alpha;
  tl2darray::reference mu;
  tl2darray::reference SiRiSr;

  tl2darray::reference alpha0;
  tl2darray::reference mu0;

  tl2darray::reference alpha1;
  tl2darray::reference mu1;

  const c_mdarray betahat;
  const c_mdarray se;
  const c_mdarray sesquare;
  const c_mdarray q;
  const c_m2darray SiRiS;

  
  tldarray::reference logodds;
  tldarray::reference sigma_beta;
  tldarray::reference sigma_beta_square;
  tl2darray::reference s;
  tl2darray::reference ssrat;

  tldarray::reference lnZ;
  tldarray::reference lnZ0;
  tldarray::reference lnZ00;
  tldarray::reference alpha_mean;
  tldarray::reference mu_mean;
  tldarray::reference pvevec;
  
  tliarray::reference iternum;
  tliarray::reference final_btnum;
  
  tldarray::reference final_max_err;
  tldarray::reference rvec;
  tldarray::reference mtp;

 
  std::vector<int> forward_range;
  std::vector<int> reverse_range;

  const size_t max_backtrack=12;
  rssr(Fit_wrap *alphas_,
       Fit_wrap *mus_,
       Fit_wrap *sirisrs_,
       const Data_wrap *datas_,
       Param_wrap *params_,
       Fit_res *results_,
       const double tolerance_, const size_t n_,
       const size_t itermax_,
       const std::vector<int> &forward_range_,
       const std::vector<int> &reverse_range_
       );

  void calc_lnZ(const blocked_range<size_t> & r)const ;
  void calculate_mtp(const blocked_range <size_t > &r)const;
  void rss_vb_iter(const blocked_range<size_t> r,bool reverse)const;
  void squarem_step(const blocked_range<size_t> & r)const;
  void squarem_backtrack(const blocked_range<size_t> & r)const;
  void calc_err(const blocked_range<size_t> &r,double &max_err)const;
  void compute_pve(const blocked_range<size_t> &r)const;
  void operator()(    const blocked_range<size_t> &r)const;
  void tls_init(const blocked_range<size_t> &r)const;
  
};





#endif
