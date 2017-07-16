norm_lfun <- function(sigu,dvec,uh,rinv,Q){
  varu <- sigu^2
  return(-1/2*(2*log(sigu)+sum(log(dvec+1/(varu)))+
                 1/(varu)%*%t(uh)%*%rinv%*%Q%*%
                 diag((varu)/(dvec*varu+1))%*%t(Q)%*%uh))
}

evdi <- function(revd){
  return(revd$vectors%*%diag(revd$values)%*%t(revd$vectors))
}

rssr_norm <- function(R,betahat,se,sigb_bounds){
  revd <- eigen(R)
  rd <- revd$values
  Q <- revd$vectors
  uh <- betahat/se
  rinv <- evdi(revd)
  ldat  <- optimise(norm_lfun,
                    interval=sigb_bounds,
                    dvec=rd,
                    uh=uh,
                    rinv=rinv,q=Q)
  return(ldat)
}






