function [pi_mean,lw] = grid_rss_varbvsr_logodds(betahat, se, SiRiS_f, sigb, log10oddsvec, alpha0,mu0)

Si=1./se(:);
SiRiS=sparse(SiRiS_f);
logoddsvec=log(10) * log10oddsvec;
[t,nl] = size(logoddsvec);
lnzl=zeros(nl,1);
for i=1:nl
%fprintf('%3f \n',logoddsvec(i));
[logw,alpha,mu,s]=rss_varbvsr(betahat,se,SiRiS,sigb,logoddsvec(i));
lnzl(i)=logw;
end
lw=normalizelogweights(lnzl);
pi=(1+power(10,-log10oddsvec)).^-1;
pi_mean=dot(pi,lw);
end