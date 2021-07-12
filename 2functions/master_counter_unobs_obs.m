function [NAMES_aux_sorted, REGION_aux_sorted, CIs,CIs_sum,post_mean_hat_sorted] = master_counter_unobs_obs(NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,policy_vars,RESULTS,confidence_level)   
%% master_counter_unobs_obs:  This function creates counterfactuals for subsets of agencies using alternative unobserved and observed covariates.  
%
% Author: José Luis Montiel Olea. Last Revised: May 26th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [NAMES_aux_sorted,...
%          REGION_aux_sorted,...
%          CIs,...
%          CIs_sum,...
%          NAIVE_counter] = master_counter_observables(NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,confidence_level)
%
% INPUT:
% NAMES_aux: Names of the agencies of interest
% Y_aux,X_aux,Z_aux: Data of the agencies of interest
% policy_vars: Vector of policy-relevant observed covariates
% RESULTS: Structure with estimation results
% confidence_level: Confidence Level

%%

K  = size(X_aux,3);  %Number of time-varying covariates

L  = size(Z_aux,2); %Number of time-invariant covariates

%% 1) Posterior Mean of alphas

addpath('../2functions')

[post_mean_hat] = posteriormean_alpha(Y_aux,X_aux,Z_aux,...
                                   RESULTS.betahat,...
                                   RESULTS.gammahat);
                               
%% 2) Counterfactuals

J_LEAs   = size(post_mean_hat,1);

COUNTERS_total = zeros(J_LEAs,J_LEAs); 
%(i,j entry: alpha of agency i, covariates of agency j)

for i_LEA = 1:J_LEAs
               
    for j_LEA = 1: J_LEAs
        
        Z_aux_policy = Z_aux(i_LEA,:);
        Z_aux_policy(1,policy_vars) = Z_aux(j_LEA, policy_vars);
        
        COUNTERS_total(i_LEA,j_LEA) ...
                 = counterfactual_function(post_mean_hat(j_LEA),...
                                        X_aux(i_LEA,:,:),...
                                        Z_aux_policy(:,:),...
                                        RESULTS.betahat(:,end),...
                                        RESULTS.gammahat(:,end))-(i_LEA==j_LEA);                                   
    end
end                                   
%% 3) Sort Counterfactuals according to Totals

TOTALS ...
                  = sum(COUNTERS_total)'; 
[~,I] ...
                  = sort(TOTALS,'descend');

NAMES_aux_sorted ...
                  = NAMES_aux(I);

COUNTERS_sorted ...
                  = COUNTERS_total(I,I); 
                           
post_mean_hat_sorted ...
                  = post_mean_hat(I);
                  
Y_aux_sorted      = Y_aux(I,:);

X_aux_sorted      = X_aux(I,:,:);

Z_aux_sorted      = Z_aux(I,:);

REGION_aux_sorted = REGION_aux(I,:);
                               
%% 4) Parametric bootstrap CIs for counters

seed        = rng('default');

Asy_Cov     =  RESULTS.asy_cov(:,:)./RESULTS.sample_sizes(end);

vector_means ...
            = [RESULTS.betahat(:,end);...
               RESULTS.gammahat(:,end)]';
                  
beta_boots  = mvnrnd(vector_means',(Asy_Cov + Asy_Cov')/2,1000);

COUNTERS_boots ...
            = zeros([size(COUNTERS_sorted),size(beta_boots,1)]);

tic;        
        
for i_boots = 1:size(beta_boots,1)
     
    beta_gamma_aux ...
            = beta_boots(i_boots,:)';
        
    beta_aux ...
            = beta_gamma_aux(1:K,:);
        
    gamma_aux ...
            = beta_gamma_aux(K+1:end,:);
        
    post_mean_hat_aux = posteriormean_alpha(Y_aux_sorted,...
                                            X_aux_sorted,...
                                            Z_aux_sorted,...
                                            beta_aux,...
                                            gamma_aux);    
    
    
    for i_LEA = 1:J_LEAs
               
        for j_LEA = 1: J_LEAs
        
        Z_aux_policy_sorted = Z_aux_sorted(i_LEA,:);
        Z_aux_policy_sorted(1, policy_vars) = Z_aux_sorted(j_LEA, policy_vars);   
         
        COUNTERS_boots(i_LEA,j_LEA,i_boots) ...
                 = counterfactual_function(post_mean_hat_aux(j_LEA),...
                                        X_aux_sorted(i_LEA,:,:),...
                                        Z_aux_policy_sorted(:,:),...
                                        beta_aux,...
                                        gamma_aux)-(i_LEA==j_LEA);                                   
        end


    end
               
end

toc;

%% 5) Quantiles for the individual effects

alpha   = 1-confidence_level;

lowerCI = quantile(COUNTERS_boots,(alpha/2),3);

upperCI = quantile(COUNTERS_boots,1-(alpha/2),3);

CIs     = cell(J_LEAs,J_LEAs);

for i_LEA = 1:J_LEAs
               
    for j_LEA = 1: J_LEAs
        
        CIs{i_LEA,j_LEA} ...
              = strcat('[',...
                        num2str(round(lowerCI(i_LEA,j_LEA),0)),...
                        ',',...
                        num2str(round(upperCI(i_LEA,j_LEA),0)),...
                        ']'); 
        
    end
    
end

%% 6) Quantiles for total effects

lowerCI_sum ...
        = reshape(quantile(sum(COUNTERS_boots,1),alpha/2,3),...
          [1,size(NAMES_aux_sorted,1)]);
      
upperCI_sum ...
        = reshape(quantile(sum(COUNTERS_boots,1),1-(alpha/2),3),...
          [1,size(NAMES_aux_sorted,1)]);  
      
CIs_sum = cell(1,J_LEAs);

for i_LEA = 1:J_LEAs 

    CIs_sum{1,i_LEA} ...
          = strcat('[',...
                        num2str(round(lowerCI_sum(1,i_LEA),0)),...
                        ',',...
                        num2str(round(upperCI_sum(1,i_LEA),0)),...
                        ']'); 
    
end

end
