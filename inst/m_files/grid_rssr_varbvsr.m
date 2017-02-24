function [lnzl] = grid_rssr_varbvsr(betahat, se, SiRiS_f, sigb, log10oddsvec, alpha0,mu0)

Si=1./se(:);
SiRiS=sparse(SiRiS_f);
logoddsvec=log(10) * log10oddsvec;
[t,nl] = size(logoddsvec);
[t,ns] = size(sigb);
lnzl=zeros(nl,ns);
for i=1:nl
    for j=1:ns
        [logw,alpha,mu,s]=rss_varbvsr(betahat,se,SiRiS,sigb(j),logoddsvec(i));
lnzl(i,j)=logw;
    end
end