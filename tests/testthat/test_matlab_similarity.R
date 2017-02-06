library(rssr)
context("Matlab concordance")

datafiles <- paste0("/media/nwknoblauch/Data/1kg/1kg_",19:22,".mat")
alpha_input <- scan("alpha0_19_22.txt",what = numeric())
mu_input <- scan("mu0_19_22.txt",what=numeric())
sigb=0.058
logodds=-2.9/log(10)
options <- list(alpha=alpha_input,mu=mu_input)
stime <- system.time(soutput <- rss_varbvsr_bigmem_squarem(datafiles,sigb=sigb,logodds=logodds,options=options))

alpha_output <- soutput$alpha
mu_output <- soutput$mu
r_output <- as.numeric(alpha_output*mu_output)

matlab_alpha_output <- scan("matlab_alpha_19_22.txt",what=numeric())
matlab_mu_output <- scan("matlab_mu_converged_19_22.txt",what=numeric())
matlab_r_output <- matlab_alpha_output*matlab_mu_output

test_that("parallel_future is somewhat comparable to MATLAB implementation",
          testthat::expect_equal(matlab_r_output,r_output,1e-5))

