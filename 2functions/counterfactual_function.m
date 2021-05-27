function [counterfactual] = counterfactual_function(alpha,X,Z,betas,gammas)
%% counterfactual_function: This function generates conuterfactuals for a given (alpha, X, Z)
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease  
%   
% Syntax: [counterfactual_function] = counterfactuals(alpha,X,Z,betas,gammas)
% 
% INPUT:
% X: 1 x T x K array of K time-varying covariates
% Z: 1 x L+1 matrix of L time-varying covariates (plus a constant)
% betas: K times 1 vector estimated coefficients
% gammas: L+1 times 1 vector of estimated coefficients

%% 1) Generate the counterfactual

K  = size(X,3); 

T  = size(X,2);
      
betahat   = betas;

gammahat  = gammas(2:end);

sum_exp_x_beta ...
          = sum(exp(reshape(X(1,:,:),[T,K])*betahat),1);

sum_exp_z_gamma ...
          = exp(Z(1,2:end)*gammahat);
    
counterfactual ...
          = alpha.* ...
            (sum_exp_z_gamma...
            *sum_exp_x_beta);    
    
end