context("prep_rss")

check_h5 <- function() {
  if (!requireNamespace("h5")) {
    skip("h5 not available")
  }
}

data("betahat")
betahat <- c(betahat)
data("se")
se <- c(se)
data("alpha_test")
data("mu_test")
data("R_shrink")
mu_test <- c(mu_test)
alpha_test <- c(alpha_test)
th5file <- system.file("example2_input.h5",package="rssr")


test_that("parameters values are the same, whether specified from a file or from a vector",{
  check_h5()

  # R_panel <-as.matrix(R_panel)
  t_SiRiS <- SiRSi(R_shrink,1/se)
#  SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
  p <- length(betahat)
  SiRiSr=(t_SiRiS%*%(alpha_test*mu_test))@x
  all_data_opts <- list(SiRiS = t_SiRiS,
                        betahat =betahat,
                        se = se,
                        alpha = alpha_test,
                        mu = mu_test,
                        tolerance= 1e-4,
                        itermax=100,
                        verbose=T,
                        lnz_tol = T)
  all_data <- prep_rss(options = all_data_opts)
  all_file_opts <- list(datafile=th5file,
                        tolerance=1e-4,
                        itermax=100,
                        verbose=T,
                        lnz_tol=T)
  all_file <- prep_rss(options=all_file_opts)
  
  expect_equal(all_file$alpha,all_data$alpha)
  expect_equal(all_file$mu,all_data$mu)
  expect_equal(all_file$SiRiS,all_data$SiRiS)
  expect_equal(all_file$betahat,all_data$betahat)
  expect_equal(all_file$se,all_data$se)
  expect_equal(all_data$sigb,all_file$sigb)
  expect_equal(all_data$logodds,all_file$logodds)
  expect_equal(all_data$SiRiSr,all_file$SiRiSr)
  
  all_file_res <- rss_varbvsr(all_file)
  all_data_res <- rss_varbvsr(all_data)
  expect_equal(all_data_res$alpha,all_file_res$alpha)
  expect_equal(all_data_res$mu,all_file_res$mu)
  expect_equal(all_data_res$SiRiSr,all_file_res$SiRiSr)
  expect_equal(all_data_res$max_err,all_file_res$max_err)
  expect_equal(all_data_res$lnZ,all_file_res$lnZ)
  expect_equal(all_data_res$iter,all_file_res$iter)
})

