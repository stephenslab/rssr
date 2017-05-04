function [betahat, se] = single_linreg(y, X)
% USEAGE: variable-by-variable regression; n obs and p vars;
% INPUT:	
%	y: n*1 vector; response; 
%	X: n*p matrix; design matrix; 
% OUTPUT: 	
%	betahat: p*1 vector; marginal effect estimate;
%	se: p*1 vector; mle std error of marginal effect;

% NB: y must be centered and X must be column centered; THIS IS IMPORTANT
	[n, p] 	 = size(X);
	SX2 	 = sum(X .* X); 
	betahat  = (X'* y)./ SX2';
	yrep 	 = repmat(y, 1, p);
	brep 	 = repmat(betahat', n, 1);
	resmat 	 = yrep - X .* brep;
	sigmahat = sqrt( sum(resmat .* resmat) / n);
	se 	 = (SX2').^(-0.5) .* sigmahat';
	betahat  = betahat(:);
	se 	 = se(:);
end