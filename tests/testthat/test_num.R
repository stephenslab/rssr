context("numerical stability")
library(Matrix)

Rlsigmoid <- function(x){
  -(log(1+exp(-x)))
}

test_that("R  C++ work identically for logsigmoid",{
  tinput <- runif(n = 10,min = 1e-30,max = 1)
  expect_equal(Rlsigmoid(tinput),logsigmoid(tinput))
})



library(rssr)
library(Matrix)

#Load the simulated data
data("betahat")
betahat <- c(betahat)
data("se")
se <- c(se)
data("alpha_test")
data("mu_test")
data("R_shrink")
#Convert the rowvector to a column vector
mu_test <- t(t(mu_test))
alpha_test <- t(t(alpha_test))
Sample_Size <- 1458
# R_panel <-as.matrix(R_panel)

#Generate SiRiS as both dense and sparse matrices 
SiRiS <- SiRSi(R_shrink,1/se)

test_that("two ways of generating siris",{
  R <- as.matrix(R_shrink)
  siris_f <- SiRSi_d(R,1/se)
  siris_f <- SiRSi_d(R,1/se)
  
})
