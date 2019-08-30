##Generate a tiny example for testing
datafiles <- paste0("/media/nwknoblauch/Data/1kg/1kg_",1:3,".mat")



pvec <- c(300,200,100)

R_list <- list()
betahat_list <- list()
se_list <- list()
alpha_list <- list()
mu_list <- list()
for(i in 1:length(datafiles)){
  R_list[[i]] <- read_ccs_h5(datafiles[i],groupname = "R")[1:pvec[i],1:pvec[i]]
  betahat_list[[i]] <- read_vec(h5filename = datafiles[i],datapath="betahat")[1:pvec[i]]
  se_list[[i]] <-  read_vec(h5filename = datafiles[i],datapath="betahat")[1:pvec[i]]
  alpha_list[[i]] <- runif(pvec[i])
  alpha_list[[i]] <- alpha_list[[i]]/sum(alpha_list[[i]])
  mu_list[[i]] <- rnorm(pvec[i])
}
sigb_test<- 0.058
logodds_test<- -2.9/log(10)

devtools::use_data(R_list,overwrite = T)
devtools::use_data(betahat_list,overwrite=T)
devtools::use_data(se_list,overwrite=T)
devtools::use_data(alpha_list,overwrite=T)
devtools::use_data(mu_list,overwrite=T)
devtools::use_data(sigb_test)
devtools::use_data(logodds_test)

R_shrink <- read_ccs_h5("~/Downloads/genotype2.mat","shrink_R")
betahat <- read_vec("~/Downloads/example2.h5","betahat")
se <- read_vec("~/Downloads/example2.h5","se")
p <- length(betahat)
alpha_test <- ralpha(p)
mu_test <- rmu(p)

devtools::use_data(R_shrink,overwrite=T)
devtools::use_data(betahat,overwrite=T)
devtools::use_data(se,overwrite=T)
devtools::use_data(alpha_test,overwrite=T)
devtools::use_data(mu_test,overwrite=T)


data("R_shrink")
data("betahat")
data("se")
data("alpha_test")
data("mu_test")

full_R <- as.matrix(R_shrink)
full_R_hf <- "/Users/nwknoblauch/Dropbox/rss_theano/example2_input.h5"
h5filename <- "/home/nwknoblauch/Dropbox/stephens/rssr/inst/example2_input.h5"
file.remove(full_R_hf)

library(h5)
hf <- h5file(full_R_hf,'a')
beta <- c(betahat)
sed <- c(se)
a_test <- c(alpha_test)
m_test <- c(mu_test)
beta_data <- createDataSet(hf,datasetname="betahat",data=beta,chunksize=NA)
se_data <- createDataSet(hf,datasetname="se",data=sed,chunksize=NA)
a_data <- createDataSet(hf,datasetname="alpha",data=a_test,chunksize=NA,compression=0L)
m_data <- createDataSet(hf,datasetname="mu",data=m_test,chunksize=NA,compression=0L)
R_data <- createDataSet(hf,datasetname="R",data=full_R,chunksize=NA,compression=0L)
h5close(hf)



hf <- h5file(h5filename,'a')
beta <- c(betahat)
sed <- c(se)
a_test <- c(alpha_test)
m_test <- c(mu_test)
beta_data <- createDataSet(hf,datasetname="betahat",data=beta,chunksize=NA)
se_data <- createDataSet(hf,datasetname="se",data=sed,chunksize=NA)
a_data <- createDataSet(hf,datasetname="alpha",data=a_test,chunksize=NA,compression=0L)
m_data <- createDataSet(hf,datasetname="mu",data=m_test,chunksize=NA,compression=0L)
h5close(hf)



data("R_panel")
data("betahat")
data("se")
R_panelm <- as.matrix(R_panel)
write.table(R_panelm,"~/Dropbox/test_rssr/R_panel.txt",sep="\t",col.names=F,row.names=F)
write.table(c(betahat),"~/Dropbox/test_rssr/betahat.txt",sep="\n",col.names=F,row.names=F)
write.table(c(se),"~/Dropbox/test_rssr/se.txt",sep="\n",col.names=F,row.names=F)

write_ccs_h5(h5filename,spmat = R_shrink,groupname = "R")

