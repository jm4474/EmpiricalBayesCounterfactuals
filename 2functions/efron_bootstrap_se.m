function [RESULTS] = efron_bootstrap_se(nboot,Y,X,Z,time_variant_variables,time_invariant_variables)
%% EFRON_BOOTSTRAP_SE: 
%   Detailed explanation goes here
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

options = optimoptions('fminunc',...
                       'Algorithm',...
                       'trust-region',...
                       'MaxIterations',...
                       10000,...
                       'MaxFunctionEvaluation',...
                       10000*K,...
                       'StepTolerance',...
                       1e-25,...
                       'FunctionTolerance',...
                       1e-20,...
                       'OptimalityTolerance',...
                       1e-5,...
                       'SpecifyObjectiveGradient',...
                        true);

problem_beta.options ...
        = options; 
            
problem_beta.solver ...
        = 'fminunc'; 
    
problem_beta.x0 ...   
        = zeros(K,1); %Initialize zeros    
    
problem_gamma.options ...
        = options;
    
problem_gamma.x0 ...   
        = zeros(L+1,1); %Initialize
    
problem_gamma.solver ...
        = 'fminunc';     

%% Estimate beta

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
    problem_gamma.objective = @(gamma) gamma_estimation(Y_aux,X_aux,Z_aux,betahat,gamma);                          

    [gammahat]= fminunc(problem_gamma); 
    
    gammas = horzcat(gammas,gammahat);

end

display(betas);
display(gammas);

% %% 5) Clean Workspace
% 
% %Results Structure 
% 
% RESULTS.betahat   = betahat;  clear betahat
% 
% RESULTS.gammahat  = gammahat;         clear gammahat
% 
% RESULTS.se_betahat ...
%                   = standard_errors(1:K);  
%               
% RESULTS.se_gammahat ...
%                   = standard_errors(K+2:K+1+L);   clear standard_errors 
% 
% RESULTS.asy_cov   = asy_cov;          clear asy_cov
% 
% RESULTS.time_varying_covs ...
%                   = COVARIATES_NAMES; clear COVARIATES_NAMES
%               
% RESULTS.time_invariant_covs ...
%                   = COVARIATES_NAMES_TIME_INVARIANT; clear COVARIATES_NAMES_TIME_INVARIANT
%               
% RESULTS.sample_sizes ...
%                   = sample_sizes; clear sample_sizes;
%               
% RESULTS.lethal_LEAs ...
%                   = lethal_LEAS; clear lethal_LEAS;
%               
% % Optimization details
%                          
% optimization_beta.objective ...
%                   = objective_fun_beta; clear objective_fun_beta
%               
% optimization_beta.exitflag ...
%                   = exitflag_beta; clear exitflag_beta
%               
% optimization_beta.fval ...
%                   = fval_beta; clear fval_beta    
%  
% optimization_beta.grad ...
%                   = grad_beta; clear grad_beta
%               
% optimization_beta.hessian ...
%                   = hessian_beta; clear hessian_beta
%               
% optimization_beta.output ...
%                   = output_beta;  clear output_beta
%               
% optimization_beta.problem ...
%                   = problem_beta; clear problem_beta
%               
%               
% optimization_gamma.objective ...
%                   = objective_fun_gamma; clear objective_fun_gamma
%               
% optimization_gamma.exitflag ...
%                   = exitflag_gamma; clear exitflag_gamma
%               
% optimization_gamma.fval ...
%                   = fval_gamma; clear fval_gamma    
%  
% optimization_gamma.grad ...
%                   = grad_gamma; clear grad_gamma
%               
% optimization_gamma.hessian ...
%                   = hessian_gamma; clear hessian_gamma
%               
% optimization_gamma.output ...
%                   = output_gamma;  clear output_gamma
%               
% optimization_gamma.problem ...
%                   = problem_gamma; clear problem_gamma  
%               
% clear options output_aux path_root ORI9 upperbound_beta;  
% 
% save('../4Output/mat/bootstrap_estimation_results.mat')
end

