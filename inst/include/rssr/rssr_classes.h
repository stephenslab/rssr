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
  double *fit_p;
  double *fit0_p;
  const bool tls;
  double *fit1_p;

  fitdtype fitl;
  fitdtype fitl0;
  fitdtype fitl1;

  Fit_wrap(double* fit_,double* fit0_,double* fit1_,const size_t p_, const size_t gridsize_);
  Fit_wrap(double* fit_,const size_t p_, const size_t gridsize_);
  Fit_wrap(double* fit_,double* fit0_, double* fit1_,const size_t p_);
  Fit_wrap(double* fit_,const size_t p_);
  void tlresize(const size_t rsize);


};


class Siris{
public:
  const size_t p;
  double* siris_p;
  Siris(double* siris_,const size_t p_);
  Siris(Eigen::MatrixXd &siris_);
};

class Data_wrap{
public:
  const size_t p;
  const double *const betahat_p;
  const double *const se_p;
  double* const sesquare_p;
  double* const q_p;

  
  Data_wrap(const double* betahat_,const double* se_,double* sesquare_,double* q_, const size_t p_);
												   
};




class Param_wrap{
public:
  const size_t gridsize;
  const size_t p;
  double *logodds_p;
  double *sigb_p;
  const bool tls;
  double *sigma_beta_square_p;
  double *s_p;
  double *ssrat_p;

  fitdtype logodds_t;
  fitdtype sigma_beta_t;
  fitdtype sigma_beta_square_t;
  fitdtype s_t;
  fitdtype ssrat_t;

  Param_wrap(double* logodds_,double* sigb_,double* sigma_beta_square_,double* s_,double* ssrat_, const Data_wrap &obj,const size_t gridsize_);
  Param_wrap(double* sigb_,double* sigma_beta_square_,double* s_,double* ssrat_, const Data_wrap &obj,const size_t gridsize_);
  Param_wrap(double* logodds_,double* sigb_, const size_t gridsize_,const size_t p_);
  Param_wrap(double* sigb_,const size_t gridsize_,const size_t p_);
  void tlresize(const blocked_range<size_t> &r,const Data_wrap &obj);
};


class Fit_res{
public:
  const size_t gridsize;
  const size_t p;
  const bool tls;
  double *lnZ_p;
  double *lnZ0_p;
  double *lnZ00_p;
  double *alpha_mean_p;
  double *mu_mean_p;
  double *pvevec_p;
  double *final_max_err_p;
  
  int    *iternum_p;
  int    *final_btnum_p;
  
  double *rvec_p;
  double *mtp_p;



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

  
  Fit_res(double* lnZ_,
	  double* lnZ0_,
	  double* lnZ00_,
	  double* alpha_mean_,
	  double* pvevec_,
	  int* iternum_,
	  double* mu_mean_,
	  double* final_max_err_,
	  double* rvec_,
	  double* mtp_,
	  int* final_btnum_,const size_t p_,const size_t gridsize_);
  Fit_res(const size_t p_,const size_t gridsize_);
  void tlresize (const size_t rsize);


};



class rssr_norm{
public:
  mutable Fit_wrap mus;
  mutable Fit_wrap sirisrs;
  mutable Fit_res results;
  const bool tls;
  const bool use_lnztol=true;
  const size_t p;
  const size_t n;
  size_t gridsize;

  const int itermax;
  const double tolerance;

  mutable Param_wrap paraams;
  const Data_wrap datas;
  const Siris siris;





  const c_mdarray betahat;
  const c_mdarray se;
  mdarray sesquare;
  mdarray q;

  
  mutable mdarray sigma_beta;
  mutable mdarray sigma_beta_square;
  mutable m2darray s;
  mutable m2darray ssrat;
  mutable mdarray lnZ;
  mutable mdarray lnZ0;
  mutable mdarray lnZ00;
  mutable mdarray alpha_mean;
  mutable mdarray mu_mean;
  mutable mdarray pvevec;
  mutable miarray iternum;
  mutable miarray final_btnum;
  mutable mdarray final_max_err;
  mutable mdarray rvec;
  mutable mdarray mtp;
  mutable m2darray SiRiSr;
  m2darray SiRiS;
  mutable m2darray mu;
  mutable m2darray mu0;
  mutable m2darray mu1;

  



  
  std::vector<int> forward_range;
  std::vector<int> reverse_range;

  const size_t max_backtrack=12;
  

  rssr_norm(const Fit_wrap &mus_,
       const Fit_wrap &sirisrs_,
       const Data_wrap &datas_,
       const Param_wrap &params_,
       const Siris &siris_,
       const Fit_res &results_,
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
  mutable Fit_wrap alphas;
  mutable Fit_wrap mus;
  mutable Fit_wrap sirisrs;
  mutable Fit_res results;
  
  const bool use_lnztol=true;
  const bool tls;
  const size_t p;
  const size_t n;


  const int itermax;
  const double tolerance;

  mutable Param_wrap paraams;
  const Data_wrap datas;
  const Siris siris;
  size_t gridsize;

  mutable m2darray alpha;
  mutable m2darray mu;
  mutable m2darray SiRiSr;

  mutable m2darray alpha0;
  mutable m2darray mu0;

  mutable m2darray alpha1;
  mutable m2darray mu1;

  const c_mdarray betahat;
  const c_mdarray se;
  const mdarray sesquare;
  const mdarray q;
  mutable mdarray logodds;
  mutable mdarray sigma_beta;
  mutable mdarray sigma_beta_square;
  mutable m2darray s;
  mutable m2darray ssrat;

  m2darray SiRiS;
  
  mutable mdarray lnZ;
  mutable mdarray lnZ0;
  mutable mdarray lnZ00;
  mutable mdarray alpha_mean;
  mutable mdarray mu_mean;
  mutable mdarray pvevec;
  
  mutable miarray iternum;
  mutable miarray final_btnum;
  
  mutable mdarray final_max_err;
  mutable mdarray rvec;
  mutable mdarray mtp;

 
  std::vector<int> forward_range;
  std::vector<int> reverse_range;

  const size_t max_backtrack=12;
  rssr(const Fit_wrap &alphas_,
       const Fit_wrap &mus_,
       const Fit_wrap &sirisrs_,
       const Data_wrap &datas_,
       const Param_wrap &params_,
       const Siris &siris_,
       const Fit_res &results_,
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
