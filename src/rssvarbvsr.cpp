#include <RcppEigen.h>
#include "rssvarbvsr.hpp"
#include "sigmoid.hpp"
#include <math.h>

// USAGE: run a single coordinate ascent of variational update to fit RSS-BVSR
// [[Rcpp::export]]
 void rss_varbvsr_update (const double betahat,
                         const double se,
                         const double sigma_beta,
                         const Eigen::ArrayXd &SiRiS_snp,
                         Eigen::ArrayXd &SiRiSr,
                         const double SiRiSr_snp,
                         const double logodds,
                         double &alpha,
                         double &mu) {
  
  double se_square = se * se;
  double sigma_beta_square = sigma_beta * sigma_beta;
  
  // Compute the variational estimate of the posterior variance.
  double sigma_square = (se_square * sigma_beta_square) / (se_square + sigma_beta_square);
  
  // Update the variational estimate of the posterior mean.
  double r = alpha * mu;
  mu = sigma_square * (betahat / se_square + r / se_square - SiRiSr_snp);
  
  // Update the variational estimate of the posterior inclusion probability.
  double SSR = mu * mu / sigma_square;
  alpha = sigmoid(logodds + 0.5 * (log(sigma_square/(sigma_beta_square)) + SSR));
  
  // Update SiRiSr = inv(S)*R*inv(S)*r
  double r_new = alpha * mu;
  SiRiSr+=(SiRiS_snp*(r_new-r));
//  add(SiRiSr, r_new - r, SiRiS_snp, p);
}
