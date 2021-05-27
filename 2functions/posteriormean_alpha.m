function [post_alpha_hat] = posteriormean_alpha(Y,X,Z,betas,gammas)                                             
%% posteriormean_alpha: This function estimates the posterior means of alpha for J agencies, given matrices Y, X, Z, betas, gammas
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [post_alpha_hat] = posteriormean_alpha(Y,X,Z,betas,gammas)
%
% INPUT:
% Y: J x T matrix of outcome variables
% X: J x T x K array of K time-varying covariates
% Z: J x L+1 matrix of L time-varying covariates (plus a constant)
% betas: K times 1 vector estimated coefficients
% gammas: L+1 times 1 vector of estimated coefficients

%% 1) Estimator of the posterior mean of alpha

K  = size(X,3); 

T  = size(Y,2);

J ...
   = size(Y,1);

post_alpha_hat ...
          = zeros(J,1);
      
sum_exp_x_beta ...
          = zeros(J,1);
      
sum_exp_z_gamma ...
          = zeros(J,1);   
      
betahat   = betas;

gammahat  = gammas(2:end);

for j_LEA = 1:J
         
    sum_exp_x_beta(j_LEA,1) ...
          = sum(exp(reshape(X(j_LEA,:,:),[T,K])*betahat),1);
      
    sum_exp_z_gamma(j_LEA,1) ...
          = exp(Z(j_LEA,2:end)*gammahat);
    
    post_alpha_hat(j_LEA,1) ...
          = (sum(Y(j_LEA,:),2)+1)./...
            (sum_exp_z_gamma(j_LEA,1)...
            *sum_exp_x_beta(j_LEA,1));    
    
end

end