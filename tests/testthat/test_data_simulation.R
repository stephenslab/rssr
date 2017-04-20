context("Data Simulation")

library(RcppOctave)
#library(rssr)
library(Matrix)
library(testthat)
#Find the location of the .m files 
mfile <- system.file("m_files/run_install.m",package="rssr")
mdir <- system.file("m_files",package="rssr")
.CallOctave('cd',mdir)
o_source("run_install.m")


test_that("Effect sizes are created in roughly the same fashion",{
  
  p <- 100000
  tpi <- 0.05
  tsigb <- 1
  
  mat_res <- .CallOctave('effectsize_maker',c(p),c(ceiling(tpi*p),0))
  Z <- rbinom(n = p,size = 1,prob = tpi)
  beta <- numeric(p)
  beta[Z==1] <- rnorm(sum(Z),mean=0,sd=tsigb)
sum(Z==1)
sum(mat_res$gamma)
  mean(mat_res$beta[mat_res$gamma==1])
  
})



