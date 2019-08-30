
library(RcppOctave)
#library(rssr)
library(Matrix)
library(testthat)
#Find the location of the .m files 
mfile <- system.file("m_files/run_install.m",package="rssr")
mdir <- system.file("m_files",package="rssr")

#change to the directory with the .m files in Octave
.CallOctave('cd',mdir)
o_source("run_install.m")

#Load the simulated data
data("betahat")
betahat <- c(betahat)
data("se")
se <- c(se)
data("alpha_test")
data("mu_test")
data("R_shrink")
#Convert the rowvector to a column vector
mu_test <- t(t(mu_test))
alpha_test <- t(t(alpha_test))
# R_panel <-as.matrix(R_panel)

#Generate SiRiS as both dense and sparse matrices 
SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- SiRSi_d(as.matrix(R_shrink),1/se)
#SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
#SiRiS <-as(SiRiS_f,"dgCMatrix")
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))

pkgdir <- "~/Dropbox/stephens/rssr"

I <- 1:p
rI <- p:1
sigb <- 1
logodds <- -3




matlab_update <- .CallOctave('wrap_squarem_adjust',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test),1e-4)
devtools::use_data(matlab_update,overwrite=T,pkg=pkgdir)





matlab_varbvsr_update_1 <- .CallOctave('wrap_rss_varbvsr_update',SiRiS_f,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,I)
devtools::use_data(matlab_varbvsr_update_1,overwrite=T,pkg=pkgdir)



rI <- p:1
sigb <- 1
logodds <- -3
matlab_varbvsr_update_2 <- .CallOctave('wrap_rss_varbvsr_update',SiRiS_f,sigb,logodds,betahat,se,alpha_test,mu_test,SiRiSr,rI)
devtools::use_data(matlab_varbvsr_update_2,overwrite=T,pkg=pkgdir)


matlab_SiRiS = .CallOctave('gen_SiRiS',as.matrix(R_shrink),se)
devtools::use_data(matlab_SiRiS,overwrite=T,pkg=pkgdir)

matlab_varbvsr_naive_1 <- .CallOctave('wrap_rss_varbvsr_naive',t(t(betahat)),t(t(se)),SiRiS_f,sigb,-4.6,t(alpha_test),t(mu_test))
devtools::use_data(matlab_varbvsr_naive_1,overwrite=T,pkg=pkgdir)


sigb <- 1  
logodds <- -3
matlab_varbvsr_squarem_1 <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test),1e-4)



matlab_varbvsr_squarem_1_up <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test),1.0)
devtools::use_data(matlab_varbvsr_squarem_1,overwrite=T,pkg=pkgdir)



sigb <- 1  
logodds <- -3
matlab_varbvsr_squarem_1 <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test),1e-4)
matlab_varbvsr_squarem_1_up <- .CallOctave('wrap_rss_varbvsr_squarem',t(t(betahat)),t(t(se)),SiRiS_f,sigb,logodds,t(alpha_test),t(mu_test),1.0)
devtools::use_data(matlab_varbvsr_squarem_1,overwrite=T,pkg=pkgdir)


sigb <- 1
log10oddsvec <- seq(-6,-1,0.5)
logoddsvec <- log10oddsvec*log(10)
matlab_grid_logodds <- .CallOctave('grid_rss_varbvsr_logodds',t(t(betahat)),t(t(se)),SiRiS_f,sigb,log10oddsvec,t(alpha_test),t(mu_test))
devtools::use_data(matlab_grid_logodds,overwrite=T,pkg=pkgdir)

log10oddsvec <- seq(-3.1,-2.1,length.out = 5)
logoddsvec <- log10oddsvec*log(10)
sigb <- seq(0.8,1.2,length.out = 5)
matlab_grid_logodds_sigb <-  .CallOctave('grid_rssr_varbvsr',t(t(betahat)),t(t(se)),SiRiS_f,sigb,log10oddsvec,t(alpha_test),t(mu_test))
devtools::use_data(matlab_grid_logodds_sigb,overwrite=T,pkg=pkgdir)
#Call with Octave and RSSR




s_test=t(t(se*se*(sigb*sigb)/(se*se+sigb)))
matlab_betavar <- .CallOctave('betavar',alpha_test,mu_test,s_test)
devtools::use_data(matlab_betavar,overwrite=T,pkg=pkgdir)

logodds <- -3
matlab_intgamma <- .CallOctave('intgamma',logodds,alpha_test)
devtools::use_data(matlab_intgamma,overwrite=T,pkg=pkgdir)

sigb <- 1
s_test=t(t(se*se*(sigb*sigb)/(se*se+sigb)))
matlab_intklbeta <- .CallOctave('intklbeta_rssbvsr',t(t(alpha_test)),t(t(mu_test)),s_test,sigb)
devtools::use_data(matlab_intklbeta,overwrite=T,pkg=pkgdir)

