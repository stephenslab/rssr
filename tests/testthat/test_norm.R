context("Normal")





test_that("Normal approximation works when beta==betahat",{
  data("betahat_mat_norm")
  data("R_panel")
  data("se_mat_norm")
  data("tparam_df_norm")
  library(dplyr)
  R_paneld <- as.matrix(R_panel)
  sigb_range <- c(min(tparam_df_norm$tsigb),max(tparam_df_norm$tsigb))
  paramdf <- data.frame(sigb=seq(sigb_range[1]/2,sigb_range[2]*2,length.out = 100))
  ng <- ncol(betahat_mat_norm)
  p <- nrow(R_panel)
  alphap=rep(1,p)
  mup <- rmu(p)
  
  tresl <- list()
  i <- 1
  siris <- SiRSi_d(R_paneld,Si = 1/se_mat_norm[,i])
  sirisr <- c(siris%*%(alphap*mup))
  tres <- grid_search_rss_varbvsr_norm(SiRiS = siris,
                                       sigma_beta = paramdf$sigb,
                                       betahat = betahat_mat_norm[,i],
                                       se = se_mat_norm[,i],talpha0 = alphap,
                                       tmu0 = mup,tSiRiSr0 = sirisr,tolerance = 1e-3,
                                       itermax = 200,verbose = F,lnz_tol = T) %>% mutate(fgeneid=i)
  
    
    
  
  fnres <- group_by(tres,fgeneid)  %>% mutate(w=normalizeLogWeights(lnZ),pve=pve/382) %>%
    summarise(mean_alpha=sum(w*alpha_mean),
              mean_pve=sum(w*pve),
              mean_sigb=sum(w*sigb)) %>% ungroup() %>% inner_join(tparam_df_norm,by="fgeneid")
  expect_equal(fnres$tsigb,fnres$mean_sigb,tolerance=0.4)
})

test_that("norm alt grid search works the same as normal norm grid search",{
  data("betahat_mat_norm")
  data("R_panel")
  data("se_mat_norm")
  data("tparam_df_norm")
  library(dplyr)
  R_paneld <- as.matrix(R_panel)
  sigb_range <- c(min(tparam_df_norm$tsigb),max(tparam_df_norm$tsigb))
  paramdf <- data.frame(sigb=seq(sigb_range[1]/2,sigb_range[2]*2,length.out = 200))
  ng <- ncol(betahat_mat_norm)
  p <- nrow(R_panel)
  alphap=rep(1,p)
  mup <- rmu(p)
  
  tresl <- list()
  i <- 1
  siris <- SiRSi_d(R_paneld,Si = 1/se_mat_norm[,i])
  sirisr <- c(siris%*%(alphap*mup))
  tres_old <- grid_search_rss_varbvsr_norm(SiRiS = siris,
                                           sigma_beta = paramdf$sigb,
                                           betahat = betahat_mat_norm[,i],
                                           se = se_mat_norm[,i],talpha0 = alphap,
                                           tmu0 = mup,tSiRiSr0 = sirisr,tolerance = 1e-3,
                                           itermax = 200,verbose = F,lnz_tol = T)
  
  
  tres_new <- grid_search_rss_varbvsr_norm_tls(SiRiS = siris,
                                               sigma_beta = paramdf$sigb,
                                               betahat = betahat_mat_norm[,i],
                                               se = se_mat_norm[,i],
                                               mu0 = mup,SiRiSr0 = sirisr,tolerance = 1e-3,
                                               itermax = 200,lnz_tol = T)
  tres_new <- tres_new[,colnames(tres_old)]
  expect_equal(tres_new,tres_old)
  
  # mcb <- microbenchmark::microbenchmark(old=grid_search_rss_varbvsr_norm(SiRiS = siris,
  #                                                                        sigma_beta = paramdf$sigb,
  #                                                                        betahat = betahat_mat_norm[,i],
  #                                                                        se = se_mat_norm[,i],talpha0 = alphap,
  #                                                                        tmu0 = mup,tSiRiSr0 = sirisr,tolerance = 1e-3,
  #                                                                        itermax = 200,verbose = F,lnz_tol = T),
  #                                         new=grid_search_rss_varbvsr_norm_tls(SiRiS = siris,
  #                                                                            sigma_beta = paramdf$sigb,
  #                                                                            betahat = betahat_mat_norm[,i],
  #                                                                            se = se_mat_norm[,i],
  #                                                                            mu0 = mup,SiRiSr0 = sirisr,tolerance = 1e-3,
  #                                                                            itermax = 200,lnz_tol = T))
  
})


test_that("norm optim works",{
  data("betahat_mat_norm")
  data("R_panel")
  data("se_mat_norm")
  data("tparam_df_norm")
  library(dplyr)
  R_paneld <- as.matrix(R_panel)
  sigb_range <- c(min(tparam_df_norm$tsigb),max(tparam_df_norm$tsigb))
  paramdf <- data.frame(sigb=seq(sigb_range[1]/2,sigb_range[2]*2,length.out = 200))
  ng <- ncol(betahat_mat_norm)
  p <- nrow(R_panel)
  mup <- rmu(p)
  
  tresl <- list()
  i <- 1
  siris <- SiRSi_d(R_paneld,Si = 1/se_mat_norm[,i])
  sirisr <- c(siris%*%(mup))
  sigbb <- c(min(paramdf$sigb),max(paramdf$sigb))
  or <- rss_varbvsr_norm_optim(SiRiS = siris,sigbb = sigbb,
                         betahat = betahat_mat_norm[,i],
                         se = se_mat_norm[,i],
                         tmu0 = mup,tSiRiSr0 = sirisr,tolerance = 1e-3,
                         itermax = 200,lnz_tol = T)
  
 
  expect_equal(tres_new,tres_old)
  
  # mcb <- microbenchmark::microbenchmark(old=grid_search_rss_varbvsr_norm(SiRiS = siris,
  #                                                                        sigma_beta = paramdf$sigb,
  #                                                                        betahat = betahat_mat_norm[,i],
  #                                                                        se = se_mat_norm[,i],talpha0 = alphap,
  #                                                                        tmu0 = mup,tSiRiSr0 = sirisr,tolerance = 1e-3,
  #                                                                        itermax = 200,verbose = F,lnz_tol = T),
  #                                         new=grid_search_rss_varbvsr_norm_tls(SiRiS = siris,
  #                                                                            sigma_beta = paramdf$sigb,
  #                                                                            betahat = betahat_mat_norm[,i],
  #                                                                            se = se_mat_norm[,i],
  #                                                                            mu0 = mup,SiRiSr0 = sirisr,tolerance = 1e-3,
  #                                                                            itermax = 200,lnz_tol = T))
  
})
