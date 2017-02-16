context("prep_rss")


data("betahat")
betahat <- c(betahat)
data("se")
se <- c(se)
data("alpha_test")
data("mu_test")
data("R_shrink")

mu_test <- t(t(mu_test))
alpha_test <- t(t(alpha_test))
# R_panel <-as.matrix(R_panel)
t_SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))
temp_h5file <- tempfile()
temp_alpha_file <- tempfile()
temp_mu_file <- tempfile()

create_h5_data(temp_h5file,R = R_shrink,betahat = betahat,se = se)
write.table(alpha_test,file = temp_alpha_file,sep="\n",col.names=F,row.names=F)
write.table(mu_test,file = temp_mu_file,sep="\n",col.names=F,row.names=F)

tadf <- tempfile()
file.create(tadf)
all_data_opts <- list(SiRiS = t_SiRiS,
                         betahat =betahat,
                         se = se,
                         alpha = c(alpha_test),
                         mu = c(mu_test),
                         tolerance= 1e-4,
                         itermax=100,
                         verbose=T,
                         lnz_tol = T)
all_data <- prep_rss(datafile=tadf,options = all_data_opts)
all_file_opts <- list(datafile=temp_h5file,
                         alphafile=temp_alpha_file,
                         mufile=temp_mu_file,
                         tolerance=1e-4,
                         itermax=100,
                         verbose=T,
                         lnz_tol=T)
all_file <- prep_rss(datafile=temp_h5file,options=all_file_opts)
test_that("parameters values are the same, whether specified from a file or from a vector",{
  expect_equal(all_file$alpha,all_data$alpha)
  expect_equal(all_file$mu,all_data$mu)
  expect_equal(all_file$SiRiS,all_data$SiRiS)
  expect_equal(all_file$betahat,all_data$betahat)
  expect_equal(all_file$se,all_data$se)
  expect_equal(all_data$sigb,all_file$sigb)
  expect_equal(all_data$logodds,all_file$logodds)
  expect_equal(all_data$SiRiSr,all_file$SiRiSr)
})
file.remove(c(temp_alpha_file,temp_mu_file,temp_h5file))
