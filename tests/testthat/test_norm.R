context("Normal")





test_that("Normal approximation works when beta==betahat",{
  data("betahat_mat_norm")
  data("R_panel")
  data("se_mat_norm")
  data("tparam_df_norm")
  R_paneld <- as.matrix(R_panel)
  sigb_range <- c(min(tparam_df_norm$tsigb),max(tparam_df_norm$tsigb))
  paramdf <- data_frame(sigb=seq(sigb_range[1]/2,sigb_range[2]*2,length.out = 100))
  ng <- ncol(betahat_mat_norm)
  alphap=ralpha(p)
  mup <- rmu(p)
  for(i in 1:ng){
    siris <- SiRSi_d(R_paneld,Si = 1/se_mat_norm[,i])
    grid_search_rss_varbvsr_norm(SiRiS = siris,sigma_beta = paramdf$sigb,)
  }
 
})
