function [RESULTS] = efron_bootstrap_se(nboot,alpha,Y,X,Z,time_variant_variables,time_invariant_variables)
%% EFRON_BOOTSTRAP: This function produces 
%
% Author: José Luis Montiel Olea. Last Revised: June 3rd, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [RESULTS] = efron_bootstrap_se(nboot,alpha,Y,X,Z,time_variant_variables,time_invariant_variables)
%
% INPUT:
% nboot: Number of boostrap samples
% alpha: Significance level for confidence intervals
% Y,X,Z: data
% time_(in)variant_variables: names of variables
%%
T  = size(Y,2);

K  = size(X,3); %Number of time-varying covariates

L  = size(Z,2); %Number of time-invariant covariates

COVARIATES_NAMES ...
   = time_variant_variables;

COVARIATES_NAMES_TIME_INVARIANT ...
   = time_invariant_variables;

%% 1) Add path for negative log-likelihood

addpath('../2functions')
%% 2) Define the optimization problem

options = optimoptions('fminunc','Algorithm','trust-region','MaxIterations',10000,'MaxFunctionEvaluation', 10000*K,...
                       'StepTolerance',...
                       1e-25,...
                       'FunctionTolerance',...
                       1e-20,...
                       'OptimalityTolerance',...
                       1e-6,...
                       'SpecifyObjectiveGradient',...
                        true, 'display', 'off');

problem_beta.options = options; 
            
problem_beta.solver = 'fminunc'; 
    
problem_beta.x0 = zeros(K,1); %Initialize zeros    
    
problem_gamma.options = options;
    
problem_gamma.x0 = zeros(L+1,1); %Initialize
    
problem_gamma.solver = 'fminunc';     

%% 3) Estimate Boostraped Parameters

tic;
Y_aux = Y; X_aux = X; Z_aux = Z;  %Pooled estimation   
    
sample_sizes = size(Y_aux,1);

betas = [];
gammas = [];

for i=1:nboot
    
    idx = randi(sample_sizes,[1,sample_sizes]);

    %Estimate Beta
    Y_aux_sample = Y_aux(idx,:);
    X_aux_sample = X_aux(idx,:);
    Z_aux_sample = Z_aux(idx,:);

    problem_beta.objective = @(beta) neg_loglikelihood(Y_aux_sample,X_aux_sample,beta);

    [betahat]= fminunc(problem_beta);

    betas = horzcat(betas,betahat);
    
    %Estimate Gamma
    problem_gamma.objective = @(gamma) gamma_estimation(Y_aux_sample,X_aux_sample,Z_aux_sample,betahat,gamma);                          

    [gammahat]= fminunc(problem_gamma); 
    
    gammas = horzcat(gammas,gammahat);
    
    display(round((i/nboot*100),2))

end

%% 4) Construct CIs
betas_CI = zeros(size(betas,1),2);
gammas_CI = zeros(size(gammas,1),2);

for i=1:size(betas,1)
    betas_CI(i,1) = quantile(betas(i,:),(alpha/2));
    betas_CI(i,2) = quantile(betas(i,:),(1-(alpha/2)));
end

for i=1:size(gammas,1)
    gammas_CI(i,1) = quantile(gammas(i,:),(alpha/2));
    gammas_CI(i,2) = quantile(gammas(i,:),(1-(alpha/2)));
end

%% 5) Calculate Standard Errors
betas_SE = zeros(size(betas,1),1);
gammas_SE = zeros(size(gammas,1),1);

for i=1:size(betas,1)
    betas_SE(i,1) = std(betas(i,:));
end

for i=1:size(gammas,1)
    gammas_SE(i,1) = std(gammas(i,:));
end

%% 6) Calculate Mean
betas_mean = zeros(size(betas,1),1);
gammas_mean = zeros(size(gammas,1),1);

for i=1:size(betas,1)
    betas_mean(i,1) = mean(betas(i,:));
end

for i=1:size(gammas,1)
    gammas_mean(i,1) = mean(gammas(i,:));
end

%% 7) Clean Workspace

RESULTS.betas   = betas; clear betas

RESULTS.gammas  = gammas; clear gammas

RESULTS.betas_CI = betas_CI; clear betas_CI
              
RESULTS.gammas_CI = gammas_CI;  clear gammas_CI 

RESULTS.betas_SE = betas_SE; clear betas_SE

RESULTS.gammas_SE = gammas_SE; clear gammas_SE

RESULTS.betas_mean = betas_mean; clear betas_mean

RESULTS.gammas_mean = gammas_mean; clear gammas_mean

RESULTS.time_variant_variables = COVARIATES_NAMES; clear COVARIATES_NAMES
              
RESULTS.time_invariant_covs = COVARIATES_NAMES_TIME_INVARIANT; clear COVARIATES_NAMES_TIME_INVARIANT
              
RESULTS.sample_sizes = sample_sizes; clear sample_sizes

save('../4Output/mat/bootstrap_estimation_results.mat')
end

