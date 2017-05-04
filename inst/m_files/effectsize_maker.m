function [beta, gamma] = effectsize_maker(p, betatype)
% USAGE: generate effect sizes for various scenarios
% INPUT:
%	p: scalar, the total number of snps analyzed
%	betatype: [num_large, num_small]
% OUTPUT:
%       beta: p by 1, true multi-snp genetic effect
%       gamma: p by 1, causal indicator for each snp (0.01 for small effect)

	num_large = betatype(1);
	num_small = betatype(2);

	beta  = zeros(p, 1);
	gamma = zeros(p, 1);

	II 	 = randperm(p); 
	I 	 = II(1:num_large);
	beta(I)  = randn(num_large, 1);
	gamma(I) = 1;
	
	if num_small ~= 0
		LL 	  = II((num_large+1):end);
		JJ 	  = LL(randperm(p-num_large)); 
		num_small = min(p-num_large, num_small);
		J 	  = JJ(1:num_small);
		beta(J)   = 0.001 * randn(num_small, 1);
		gamma(J)  = 0.001;
	end

	if (sum(gamma == 1) ~= num_large) || (sum(gamma == 0.001) ~= num_small)
		error('check the generation of effect sizes')
	end
end
