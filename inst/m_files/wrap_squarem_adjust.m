function [mtp, alpha,alpha0,alpha1,alpha2, mu,mu0,mu1,mu2,SiRiSr] = wrap_squarem_adjust(betahat, se, SiRiS, sigb, logodds, alpha0,mu0,tolerance)


options=struct('alpha',alpha0,'mu',mu0,'tolerance',tolerance);
% Get the number of variables (p).
p = length(betahat);

% SiRiS must be a sparse matrix.
if ~issparse(SiRiS)
    SiRiS = sparse(double(SiRiS));
end

%  if ~exist('options', 'var')
%  options = [];
%end

if isfield(options,'tolerance')
    tolerance = double(options.tolerance);
else
    tolerance = 1e-4; 
end

% Set initial estimates of variational parameters.
if isfield(options,'alpha')
    alpha = double(options.alpha(:));
else
    alpha = rand(p,1);
    alpha = alpha / sum(alpha);
end
if isfield(options,'mu')
    mu = double(options.mu(:));
else
    mu = randn(p,1);
end
if length(alpha) ~= p || length(mu) ~= p
    error('options.alpha and options.mu must be vectors of the same length');
end

% Determine whether to update the prior SD of the additive effects.
if isfield(options,'update_sigb')
    update_sigb = options.update_sigb;
else
    update_sigb = false;
end

% Determine whether to display the algorithm's progress.
if isfield(options,'verbose')
    verbose = options.verbose;
else
    verbose = true;
end

% Determine whether to modify the step length in SQUAREM (step 6 in Table 1).
if isfield(options,'modify_step')
    modify_step = options.modify_step;
else
    modify_step = true;
end

clear options;

% Compute a few useful quantities for the main loop.
SiRiSr = full(SiRiS * (alpha .* mu));
q 	 = betahat ./ (se .^2);

% Calculate the variance of the coefficients.
se_square 	= se .* se;
sigb_square 	= sigb * sigb;
s 		= (se_square .* sigb_square) ./ (se_square + sigb_square);

% Initialize the fields of the structure info.
iter   = 0;
loglik = [];

% Calculate the variational lower bound based on the initial values.
r   = alpha .* mu;
lnZ = q'*r - 0.5*r'*SiRiSr - 0.5*(1./se_square)'*betavar(alpha, mu, s);
lnZ = lnZ + intgamma(logodds, alpha) + intklbeta_rssbvsr(alpha, mu, s, sigb_square);
%fprintf('Calculate the variational lower bound based on the initial values: %+13.6e ...\n', lnZ);

loglik = [loglik; lnZ]; %#ok<AGROW>

if verbose
    %fprintf('       variational    max. incl max.       \n');
    %fprintf('iter   lower bound  change vars E[b] sigma2\n');
end

% Repeat until convergence criterion is met.

% Go to the next iteration.
iter = iter + 1;

% Save the current variational parameters and lower bound.
alpha0  = alpha; 
mu0     = mu;
lnZ0    = lnZ;
params0 = [alpha; alpha .* mu];

% Run a forward or backward pass of the coordinate ascent updates.
if mod(iter,2)
    I = (1:p);
else
    I = (p:-1:1);
end
% %fprintf('       First Update       \n');

[alpha1, mu1, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha0, mu0, SiRiSr, I); 
%      %fprintf('       Second Update   \n');

% Run the second fix-point mapping step (line 2 of Table 1).
[alpha2, mu2, SiRiSr] = rss_varbvsr_update(SiRiS, sigb, logodds, betahat, se, alpha1, mu1, SiRiSr, I);

% Compute the step length (line 3-5 of Table 1).
alpha_r 	= alpha1 - alpha0;
mu_r 	= mu1 - mu0;
alpha_v 	= (alpha2 - alpha1) - alpha_r;
mu_v 	= (mu2 - mu1) - mu_r;
mtp 	= - sqrt(norm(alpha_r)^2+norm(mu_r)^2) / ...
          sqrt(norm(alpha_v)^2+norm(mu_v)^2);

  
   
				 % set mtp = -1
    % scenario 2: no need to modify the step length (line 7 of Table 1)  
    fprintf('Squarem_adjust!:       %d   \n',iter);% i.e. mtp < -1
    alpha  = alpha0 - 2*mtp*alpha_r + (mtp^2)*alpha_v;
    mu     = mu0 - 2*mtp*mu_r + (mtp^2)*mu_v;
    %    SiRiSr3 = full(SiRiS * (alpha3 .* mu3));




end
