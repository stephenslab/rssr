## ---- eval=FALSE---------------------------------------------------------
#  devtools::install_github("stephenslab/rssr",ref="v0.1.0-alpha")

## ---- message=FALSE------------------------------------------------------
suppressPackageStartupMessages(library(rssr))
suppressPackageStartupMessages(library(Matrix))

# marginal effect estimates 
data("betahat")
betahat <- c(betahat)

# standard errors of marginal effect estimates
data("se")
se <- c(se)

# shrinkage estimate of correlation matrix
# note that this matrix is sparse: look at `is(R)`
data("R_shrink")
R <- R_shrink

# initial values of variational parameters 
data("alpha_test")
data("mu_test")

# convert row vectors to column vectors
mu_start <- t(t(mu_test))
alpha_start <- t(t(alpha_test))

## ------------------------------------------------------------------------
# compute the matrix SiRiS
SiRiS_f <- as.matrix(rssr::SiRSi(R, 1/se))
# note that SiRiS can be also computed as follows:
# SiRiS_f <- diag(1/se) %*% as.matrix(R) %*% diag(1/se)

# convert SiRiS to a sparse matrix
SiRiS <- as(SiRiS_f,"dgCMatrix")

# compute the initial value of SiRiSr
SiRiSr_start <- c(SiRiS_f %*% (alpha_start * mu_start))

## ------------------------------------------------------------------------
r_results <- rssr::rss_varbvsr_naive(SiRiS, sigma_beta=1, logodds=-4.6, betahat, se, alpha_start, mu_start, SiRiSr_start,itermax=100,verbose = F,lnz_tol = F, tolerance=1e-4)

# extract certain quantities for further comparisons
r_iter <- r_results$iter
r_maxerr <- r_results$max_err
r_lnZ <- r_results$lnZ
r_alpha <- r_results$alpha
r_mu <- r_results$mu

## ------------------------------------------------------------------------
matfile_name <- system.file("example2_vb_results.mat",package="rssr")
m_results <- R.matlab::readMat(matfile_name)

# number of iterations
m_iter <- c(m_results$info[[1]])

# maximum relative difference between the parameters at the last two iterations
m_maxerr <- c(m_results$info[[2]])

# variational lower bound at the last iteration
m_lnZ <- c(m_results$logw)

# variational parameter estimates
m_alpha <- c(m_results$alpha)
m_mu <- c(m_results$mu)

## ---- echo=FALSE---------------------------------------------------------
tab.mat <- matrix(data=0, nrow=3, ncol=2)
tab.mat[, 1] <- c(r_iter, log10(r_maxerr), r_lnZ)
tab.mat[, 2] <- c(m_iter, log10(m_maxerr), m_lnZ)

rownames(tab.mat) <- c("`iter`", "`log10(maxerr)`", "`lnZ`")
colnames(tab.mat) <- c("`rss_varbvsr_naive`", "`rss_varbvsr.m`")

suppressPackageStartupMessages(library(knitr))
knitr::kable(tab.mat)

## ---- echo=FALSE, fig.width=8, fig.height=4------------------------------
par(mfrow=c(1,2))

plot(r_mu, m_mu, pch=20, xlab="rss_varbvsr_naive", ylab="rss_varbvsr.m", main="mu")
abline(0,1,col="red")

plot(r_alpha, m_alpha, pch=20, xlab="rss_varbvsr_naive", ylab="rss_varbvsr.m", main="alpha")
abline(0,1,col="red")

## ---- echo=FALSE---------------------------------------------------------
sessionInfo()

