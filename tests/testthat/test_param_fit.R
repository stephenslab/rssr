context("Fitting parameters")

data(betahat)

data(se)

data(alpha_test)
data(mu_test)
data(R_shrink)

betahat <- c(betahat)
se <- c(se)
mu_test <- c(mu_test)
alpha_test <- c(alpha_test)
t_SiRiS <- SiRSi(R_shrink,1/se)
SiRiS_f <- as.matrix(SiRSi(R_shrink,1/se))
p <- length(betahat)
SiRiSr=c(SiRiS_f%*%(alpha_test*mu_test))


sigb <- 1
log10oddsvec <- seq(-6,-1,0.5)
logodds <- log10oddsvec*log(10)
# input_h5files <- dir("/media/nwknoblauch/Data/GTEx/1kg_LD",full.names = T)[1]
# output_h5files <- gsub("1kg_LD","1kg_IBD",input_h5files)
# file.remove(output_h5files)

SiRiS  <-as(SiRiS_f,"dgCMatrix")
toutfile <- tempfile()
options <- list(plan=list(engine="MC",resources=list(nodes=2)),
                betahat=betahat,
                se=se,
                logodds=logodds,
                toFile=toutfile,
                alpha=alpha_test,
                mu=mu_test,
                SiRiS=SiRiS,
                sigb=sigb,
                datafile=tempfile(),
                tolerance=1e-4)
datafiles <- tempfile()
mmat <- rss_varbvsr_parallel_future(options=options)
lnzmat <- gen_lnzmat(mmat[[1]],logodds ,sigb)
test_that("We can estimate a known value of pi by grid search",
          expect_equal(c(marg_pi(log10oddsvec,lnzmat)),(10/p),tolerance=0.01)
)
sigb <- seq(0.5,1.5,0.1)
log10oddsvec <- log((10/p)/(1-(10/p)),base = 10)
logodds <- log10oddsvec*log(10)
options <- list(plan=list(engine="MC",resources=list(nodes=2)),
                betahat=betahat,
                se=se,
                logodds=logodds,
                toFile=toutfile,
                alpha=alpha_test,
                mu=mu_test,
                SiRiS=SiRiS,
                sigb=sigb,
                datafile=tempfile(),
                tolerance=1e-4)
datafiles <- tempfile()
mmat <- rss_varbvsr_parallel_future(options=options)
lnzmat <- gen_lnzmat(mmat[[1]],logodds ,sigb)
gsigb <- marg_param(c(lnzmat),param = sigb)
test_that("We can estimate a value of sigb by grid search",
          expect_equal(gsigb,1,tolerance=0.1))






