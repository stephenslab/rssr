// Author: Peter Carbonetto (ported to RcppEigen by Nicholas Knoblauch)
// Source: https://github.com/pcarbo/varbvs/blob/master/R/varbvs/src/sigmoid.h
#ifndef INCLUDE_SIGMOID
#define INCLUDE_SIGMOID
#include <RcppEigen.h>
#include "rssr_types.h"


// Function declarations.
// -----------------------------------------------------------------
// Compute log(1 + exp(x)) in a numerically stable manner.
// Eigen::ArrayXd logpexp (Eigen::ArrayXd x);
// 
// 
// // Return the sigmoid function at x.
// Eigen::ArrayXd sigmoid (Eigen::ArrayXd x);
// 
// Eigen::ArrayXd logsigmoid (Eigen::ArrayXd x);
// Return the logarithm of the sigmoid function at x. Computation is
// performed in a numerically stable manner.


inline double logpexp (const double x) {
  return (x >= 16) * x + (x < 16)  * log(1 + exp(x));
}

inline double sigmoid (const double x) {
  return 1/(1 + exp(-x));
}

inline Eigen::ArrayXd sigmoid (const c_arrayxd_internal x) {
  return 1/(1 + (-x).exp());
}

inline double logsigmoid (const double x) {
  return -logpexp(-x);
}

inline Eigen::ArrayXd logsigmoid(const c_arrayxd_internal x){
  return -(-x).unaryExpr(&logpexp);
}




//inline double logpexp (double x);


// Return the sigmoid function at x.
//inline double sigmoid (double x);


// Return the logarithm of the sigmoid function at x. Computation is
// performed in a numerically stable manner.





#endif
