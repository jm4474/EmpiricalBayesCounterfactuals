function [RESULTS] = estimation_beta_gamma(Y,X,Z,time_variant_variables,time_invariant_variables,dummy_region_aux,REGION,NAMES)
%% estimation_beta_gamma: This function estimates the parameters beta and gamma
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [RESULTS] = estimation_beta_gamma(Y,X,Z,...
%                                       time_variant_variables,...
%                                       time_invariant_variables,...
%                                       dummy_region,REGION,NAMES)
%
% INPUT:
% Y,X,Z: data
% time_(in)variant_variables: names of variables
% dummy_region: false if off, true if on. 
% REGION: region of each agency
% NAMES: names of agencies

%%
T  = size(Y,2);

K  = size(X,3); %Number of covariates

L  = size(Z,2); %Time-invariant covariates

COVARIATES_NAMES ...
   = time_variant_variables;

COVARIATES_NAMES_TIME_INVARIANT ...
   = time_invariant_variables;

%% 1) Add path for the negative log-likelihood function

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
                    
%,...
                      %  'HessianFcn',...
                      %  'objective'                    
                    
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
    
%% 3) Create Regional Dummies

if dummy_region_aux == true
   
[Z] = create_dummy_region(REGION,NAMES,Z);

COVARIATES_NAMES_TIME_INVARIANT(1,L+1) ={'dummy_east'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+2) ={'dummy_MW'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+3) ={'dummy_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+4) ={'OFF/POP_interact_dummy_east'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+5) ={'OFF/POP_interact_dummy_MW'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+6) ={'OFF/POP_interact_dummy_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+7) ={'POV_interact_dummy_east'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+8) ={'POV_interact_dummy_MW'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+9) ={'POV_interact_dummy_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+10) ={'BLACK_interact_dummy_east'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+11) ={'BLACK_interact_dummy_MW'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+12) ={'BLACK_interact_dummy_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+13) ={'HISPANIC_interact_dummy_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+14) ={'HISPANIC_interact_dummy_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+15) ={'HISPANIC_interact_dummy_S'};

L = size(Z,2);

problem_gamma.x0 ...   
        = zeros(L+1,1); %Initialize
    
else

end


%{
REGION ...
    = cellstr(string(REGION));

NAMES ...
    = cellstr(string(NAMES));

REGION_UNIQUE ...
    = unique(REGION);

Number_regions ...
    = size(REGION_UNIQUE,1);

dummy_region ...
    = zeros(size(NAMES,1),Number_regions-1);

for i_region = 1:Number_regions-1
   
    dummy_region(:,i_region) ...
             = strcmp(REGION_UNIQUE(i_region),REGION);
         
end

%Dummy for time-invariant

dummy_region_z_on = dummy_region_aux;

if dummy_region_z_on == true

extra = dummy_region;

Z = [Z,...     
     dummy_region(:,1),...     
     dummy_region(:,2),...
     dummy_region(:,3)];
 
%Z(:,2).*dummy_region(:,1),... 
%Z(:,2).*dummy_region(:,2),...
%Z(:,2).*dummy_region(:,3),...

%COVARIATES_NAMES_TIME_INVARIANT(1,L+1) ={'OFF/POP_interact_dummy_east'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+1) ={'dummy_east'};

%COVARIATES_NAMES_TIME_INVARIANT(1,L+3) ={'OFF/POP_interact_MW'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+2) ={'dummy_MW'};

%COVARIATES_NAMES_TIME_INVARIANT(1,L+5) ={'OFF/POP_interact_S'};

COVARIATES_NAMES_TIME_INVARIANT(1,L+3) ={'dummy_S'};

L_without_dummies ...
   = L; 

L = size(Z,2);

problem_gamma.x0 ...   
        = zeros(L+1,1); %Initialize
    
else

L_without_dummies ...
   = L;     
    
end   

%}

%% 4) Initialize variables for estimation

standard_errors ...
    = zeros(K+L+1,1);

asy_cov ...
    = zeros(K+L+1, K+L+1, 1 );

%Optimization output

output_beta ...
    = struct([]);

output_gamma ...
    = struct([]);

%% 5.1) Estimate beta
% (pooled) 

tic;

Y_aux = Y; X_aux = X; Z_aux = Z;  %Pooled estimation   
    
sample_sizes ...
      = size(Y_aux,1);

lethal_LEAS ... 
      = sum(sum(Y_aux,2)>0);

% Estimate beta  
objective_fun_beta ...
    = @(beta) neg_loglikelihood(Y_aux,X_aux,beta);

problem_beta.objective ...
        = objective_fun_beta;

upperbound_beta ...
        = mean(sum(Y_aux,2)).*(log(T)/T);  

[betahat,...
 fval_beta,...
 exitflag_beta,...
 output_aux,...
 grad_beta,...
 hessian_beta] ...
        = fminunc(problem_beta); 
    
output_beta = [output_beta, output_aux];    
 
clear output_aux

%% 5.2) Estimate gamma

% estimate alpha using betahat

objective_fun_gamma ...
    = @(gamma) gamma_estimation(Y_aux,X_aux,Z_aux,...
                               betahat,gamma);                          
problem_gamma.objective ...
        = objective_fun_gamma;
    
[gammahat,...
 fval_gamma,...
 exitflag_gamma,...
 output_aux,...
 grad_gamma,...
 hessian_gamma] ...
        = fminunc(problem_gamma);     

output_gamma = [output_gamma, output_aux];   

[standard_errors,asy_cov] ...
             = std_errors(Y_aux,X_aux,Z_aux,...
                          betahat,...
                          gammahat,...
                          full(hessian_beta),...
                          full(hessian_gamma));

clear X_aux Y_aux Z_aux

toc; 

%% 5) Clean Workspace

%Results Structure 

RESULTS.betahat   = betahat;  clear betahat

RESULTS.gammahat  = gammahat;         clear gammahat

RESULTS.se_betahat ...
                  = standard_errors(1:K);  
              
RESULTS.se_gammahat ...
                  = standard_errors(K+2:K+1+L);   clear standard_errors 

RESULTS.asy_cov   = asy_cov;          clear asy_cov

RESULTS.time_varying_covs ...
                  = COVARIATES_NAMES; clear COVARIATES_NAMES
              
RESULTS.time_invariant_covs ...
                  = COVARIATES_NAMES_TIME_INVARIANT; clear COVARIATES_NAMES_TIME_INVARIANT
              
RESULTS.sample_sizes ...
                  = sample_sizes; clear sample_sizes;
              
RESULTS.lethal_LEAs ...
                  = lethal_LEAS; clear lethal_LEAS;
              
% Optimization details
                         
optimization_beta.objective ...
                  = objective_fun_beta; clear objective_fun_beta
              
optimization_beta.exitflag ...
                  = exitflag_beta; clear exitflag_beta
              
optimization_beta.fval ...
                  = fval_beta; clear fval_beta    
 
optimization_beta.grad ...
                  = grad_beta; clear grad_beta
              
optimization_beta.hessian ...
                  = hessian_beta; clear hessian_beta
              
optimization_beta.output ...
                  = output_beta;  clear output_beta
              
optimization_beta.problem ...
                  = problem_beta; clear problem_beta
              
              
optimization_gamma.objective ...
                  = objective_fun_gamma; clear objective_fun_gamma
              
optimization_gamma.exitflag ...
                  = exitflag_gamma; clear exitflag_gamma
              
optimization_gamma.fval ...
                  = fval_gamma; clear fval_gamma    
 
optimization_gamma.grad ...
                  = grad_gamma; clear grad_gamma
              
optimization_gamma.hessian ...
                  = hessian_gamma; clear hessian_gamma
              
optimization_gamma.output ...
                  = output_gamma;  clear output_gamma
              
optimization_gamma.problem ...
                  = problem_gamma; clear problem_gamma  
              
clear options output_aux path_root ORI9 upperbound_beta;  

save('../4Output/mat/estimation_results.mat')

end