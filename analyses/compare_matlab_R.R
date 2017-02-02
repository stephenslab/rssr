library(testthat)
library(RSSR)
library(Matrix)
library(h5)

alpha_input <- scan("~/Downloads/alpha0.txt",what=numeric())
mu_input <- scan("~/Downloads/mu0.txt",what=numeric())
logodds_input <- scan("~/Downloads/logodds.txt",what=numeric())
betahat_input <- scan("~/Downloads/betahat.txt",what=numeric())
se_input <- scan("~/Downloads/se.txt",what=numeric())
SiRiSr_input <- scan("~/Downloads/SiRiSr.txt",what=numeric())
sigbeta_input <- scan("~/Downloads/sigma_beta.txt",what=numeric())
SiRiS_input <- gen_SiRSi("/media/nwknoblauch/Data/1kg/1kg_19.mat")


alpha_update_matlab <- scan("~/Dropbox/RSSR/analyses/alpha_update.txt",what=numeric())
mu_update_matlab <- scan("~/Dropbox/RSSR/analyses/mu_update.txt",what=numeric())
SiRiSr_update_matlab <- scan("~/Dropbox/RSSR/analyses/SiRiSr_update.txt",what=numeric())



datafiles <- paste0("/media/nwknoblauch/Data/1kg/1kg_",19:22,".mat")
sigb=0.058
logodds=-2.9/log(10)
alphavec <- scan("~/Dropbox/RSSR/analyses/alpha_input.txt",what=numeric())
muvec <- scan("~/Dropbox/RSSR/analyses/mu_input.txt",what=numeric())
options <- list(alpha=alphavec,mu=muvec)
output <- rss_varbvsr_bigmem_squarem(datafiles,sigb,logodds,options)
str(output)



update_matlf <- "~/Dropbox/RSSR/analyses/fvarbvsr_19_22.mat"
update_matl <- h5file(update_matlf,"r")
alpha_update_matlab <- c(update_matl["alpha"][])
lnZ_update_matlab <- c(update_matl["lnZ"][])




alpha_update_R <- output$alpha
lnZ_update_R <- output$lnZ












alpha_input <- output$alpha0
write.table(alpha_input,"~/Dropbox/RSSR/analyses/alpha_input.txt",sep="\n",col.names=F,row.names=F)
mu_input <- output$mu0
write.table(mu_input,"~/Dropbox/RSSR/analyses/mu_input.txt",sep="\n",col.names=F,row.names=F)





                          

system.time(ret <- rss_varbvsr(SiRiS  =SiRiS_input,
                               sigma_beta = sigbeta_input,
                               logodds = logodds_input,
                               betahat = betahat_input,
                               se = se_input,
                               alpha = alpha_input,
                               mu = mu_input,
                               SiRiSr = SiRiSr_input))

alpha_update_R <- ret[,1]
mu_update_R <- ret[,2]
SiRiSr_update_R <- ret[,3]

expect_equal(alpha_update_R,alpha_update_matlab,tolerance=1e-6)
expect_equal(mu_update_matlab,mu_update_R,tolerance=1e-6)
expect_equal(SiRiSr_update_matlab,SiRiSr_update_R)
