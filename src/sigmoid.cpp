#include "sigmoid.hpp"
#include <math.h>
#include <RcppEigen.h>
// Function definitions.
// -----------------------------------------------------------------
// Computes log(1 + exp(x)).
/*
Eigen::ArrayXd logpexp (Eigen::ArrayXd x) {
  return (x >= 16) * x + (x < 16)  * log(1 + x.exp());
}

// Return the sigmoid function at x.
Eigen::ArrayXd sigmoid (Eigen::ArrayXd x) {
  return 1/(1 + exp(-x));
}


// Return the logarithm of the sigmoid function at x.
Eigen::ArrayXd logsigmoid (Eigen::ArrayXd x) {
  return -logpexp(-x);
}
*/
// double logpexp (double x) {
//   return (x >= 16) * x + (x < 16)  * log(1 + exp(x));
// }
// 
// 
// 
// // Return the sigmoid function at x.
// double sigmoid (double x) {
//   return 1/(1 + exp(-x));
// }
// 
// 
// // Return the logarithm of the sigmoid function at x.
// double logsigmoid (double x) {
//   return -logpexp(-x);
// }




