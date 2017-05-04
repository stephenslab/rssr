context("numerical stability")
library(Matrix)

Rlsigmoid <- function(x){
  -(log(1+exp(-x)))
}

test_that("R  C++ work identically for logsigmoid",{
  tinput <- runif(n = 10,min = 1e-30,max = 1)
  expect_equal(Rlsigmoid(tinput),logsigmoid(tinput))
})
