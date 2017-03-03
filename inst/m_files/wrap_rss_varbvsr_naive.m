function [lnZ, alpha, mu, s, info] = wrap_rss_varbvsr_naive(betahat, se, SiRiS, sigb, logodds, alpha0,mu0)
options=struct('alpha',alpha0,'mu',mu0);
[lnZ, alpha, mu, s, info] = rss_varbvsr(betahat, se, SiRiS, sigb, logodds, options);
end;
