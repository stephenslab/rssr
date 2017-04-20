context("numerical stability")


Rlsigmoid <- function(x){
  -(log(1+exp(-x)))
}

test_that("R  C++ work identically for logsigmoid",{
  tinput <- runif(n = 10,min = 1e-30,max = 1)
  expect_equal(Rlsigmoid(tinput),logsigmoid(tinput))
})


Rsigmoid <- function(x){
  1/(1+exp(-x))
}
test_that("R  C++ work identically for sigmoid",{
  tinput <- runif(n = 100000,min = 1e-30,max = 1)
  # expect_equal(sigmoid(tinput),Rsigmoid(tinput))
  # expect_equal(sigmoid_alt(tinput),Rsigmoid(tinput))
  
  res <- microbenchmark::microbenchmark(alt=sigmoid_alt(tinput),old=sigmoid(tinput))
  
})
