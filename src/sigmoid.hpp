// Author: Peter Carbonetto (ported to RcppEigen by Nicholas Knoblauch)
// Source: https://github.com/pcarbo/varbvs/blob/master/R/varbvs/src/sigmoid.h
#ifndef INCLUDE_SIGMOID
#define INCLUDE_SIGMOID
#include <RcppEigen.h>

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


inline double logpexp (double x) {
  return (x >= 16) * x + (x < 16)  * log(1 + exp(x));
}

inline double sigmoid (double x) {
  return 1/(1 + exp(-x));
}

inline double logsigmoid (double x) {
  return -logpexp(-x);
}

inline Eigen::ArrayXd logsigmoid(Eigen::ArrayXd x){
  return -(-x).unaryExpr(&logpexp);
}




//inline double logpexp (double x);


// Return the sigmoid function at x.
//inline double sigmoid (double x);


// Return the logarithm of the sigmoid function at x. Computation is
// performed in a numerically stable manner.





#endif
