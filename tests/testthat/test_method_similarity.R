library(rssr)
context("Method concordance")
data("R_list")
data("betahat_list")
data("se_list")
data("alpha_list")
data("mu_list")
data("sigb_test")
data("logodds_test")

temp_datafiles <- replicate(length(R_list),tempfile())
for(i in 1:length(temp_datafiles)){
  create_h5_data(temp_datafiles[i],R_list[[i]],betahat_list[[i]],se_list[[i]],chr=i)
}

options <- list(alpha=unlist(alpha_list),mu=unlist(mu_list))
s_imp <- rss_varbvsr_bigmem_squarem(datafiles = temp_datafiles,sigb = sigb_test,logodds = logodds_test,options = options)
p_imp <- rss_varbvsr_parallel_future(datafiles = temp_datafiles,sigb = sigb_test,logodds = logodds_test,options = options)
first_imp <- rss_varbvsr_future()
