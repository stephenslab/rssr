function sigma = usepve(X, beta, n, pve)
% USAGE: find residual SD to obtain given pve under model Y=XB+E
% INPUT:
%	X: n by p genotype matrix
%	beta: p by 1 genetic effect
%	n: sample size
%	pve: prespecified pve value
% OUTPUT:
%	sigma: residual SD

	part_1 = (1-pve) / pve;
	xb     = X * beta;  
	part_2 = dot(xb, xb) / n;
	sigma  = sqrt(part_1 * part_2);
end