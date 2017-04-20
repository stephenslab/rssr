function [y, X, beta, gamma, Nsnp, sigma] = XY_maker(X, betatype, pve, myseed)
% USAGE: generate continuous phenotype under additive model
% INPUT:
%	X: n by p genotype matrix
%	betatype: [num_large, num_small]
%	pve: scalar, user-defined pve
%	myseed: integer, random seed used in data generation
% OUTPUT:
%	y: n by 1 centered trait vector
%	X: n by p column-centered genotype matrix
%	beta: p by 1, true multi-snp genetic effect
%	gamma: p by 1, causal indicator for each snp
%	Nsnp: p by 1, sample size for each snp
%	sigma: residual SD to obtain the given pve under model Y=XB+E

	rng(myseed, 'twister');
	[n, p] 		= size(X);
	Nsnp 		= n * ones(p, 1);	
	[beta, gamma] 	= effectsize_maker(p, betatype); 
        X 		= X - repmat(mean(X),n,1); 		% center columns of genotypes
	sigma 		= usepve(X, beta, n, pve);		% decide residual sd based on pve
	y 		= X * beta + sigma * randn(n,1);		
	y 		= y - mean(y); 				% center trait
end