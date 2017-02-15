function [alpha1,mu1,SiRiSr] = wrap_rss_varbvsr_update(SiRiS_f,sigb,logodds,betahat,se,alpha0,mu0,SiRiSr,p)
  SiRiS=sparse(SiRiS_f);
  I=1:p;
  [alpha1, mu1, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha0, mu0, SiRiSr, I); 
end
