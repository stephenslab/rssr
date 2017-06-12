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
  for(i in 1:ng){
    siris <- SiRSi_d(R_paneld,Si = 1/se_mat_norm[,i])
    sirisr <- c(siris%*%(alphap*mup))
    tresl[[i]] <- grid_search_rss_varbvsr_norm(SiRiS = siris,
                                               sigma_beta = paramdf$sigb,
                                               betahat = betahat_mat_norm[,i],
                                               se = se_mat_norm[,i],talpha0 = alphap,
                                               tmu0 = mup,tSiRiSr0 = sirisr,tolerance = 1e-3,
                                               itermax = 200,verbose = F,lnz_tol = T) %>% mutate(fgeneid=i)
    
  }
  tres <- bind_rows(tresl)
  fnres <- group_by(tres,fgeneid)  %>% mutate(w=normalizeLogWeights(lnZ),pve=pve/382) %>%
    summarise(mean_alpha=sum(w*alpha_mean),
              mean_pve=sum(w*pve),
              mean_sigb=sum(w*sigb)) %>% ungroup() %>% inner_join(tparam_df_norm,by="fgeneid")
  expect_equal(fnres$tsigb,fnres$mean_sigb,tolerance=0.4)
 
})
