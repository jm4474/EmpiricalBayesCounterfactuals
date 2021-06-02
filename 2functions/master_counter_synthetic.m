function [CIs, CIs_sum] =  master_counter_synthetic(Y_aux,X_aux,Z_aux,RESULTS,confidence_level)
%% master_counter_synthetic: This function creates counterfactuals for subsets of agencies.
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [CIs, CIs_sum] =  master_counter_synthetic(Y_aux,X_aux,Z_aux,RESULTS,confidence_level)
%
% INPUT:
% Y_aux, X_aux, Z_aux: Data of the agencies of interest
% RESULTS: Structure with estimation results
% confidence_level: Confidence Level

%%
K  = size(X_aux,3);  %Number of time-varying covariates

J_LEAs ...
   = 4; 

rng('default');

Asy_Cov     =  RESULTS.asy_cov(:,:)./RESULTS.sample_sizes(end);

vector_means ...
            = [RESULTS.betahat(:,end);...
               RESULTS.gammahat(:,end)]';
                  
beta_boots  = mvnrnd(vector_means',(Asy_Cov + Asy_Cov')/2,1000);

COUNTERS_boots_1 ...
            = zeros([J_LEAs,J_LEAs,size(beta_boots,1)]);
        
COUNTERS_boots_2 ...
            = zeros([J_LEAs,J_LEAs,size(beta_boots,1)]); 
        
COUNTERS_boots ...
            = zeros([J_LEAs,J_LEAs,size(beta_boots,1)]);         
                
Aux_dummy = [1,0,0; 0,1,0; 0,0,1; 0,0,0];        

tic;        
        
for i_boots = 1:size(beta_boots,1)
     
    beta_gamma_aux ...
            = beta_boots(i_boots,:)';
        
    beta_aux ...
            = beta_gamma_aux(1:K,:);
        
    gamma_aux ...
            = beta_gamma_aux(K+1:end,:);
            
    for i_LEA = 1:J_LEAs
               
        for j_LEA = 1: J_LEAs
        
        COUNTERS_boots_1(i_LEA,j_LEA,i_boots) ...
                 = counterfactual_function(1,...
                                        X_aux(i_LEA,:,:),...
                                        [1,Z_aux(i_LEA,:),...
                                        Aux_dummy(j_LEA,:),...
                                        Z_aux(i_LEA,2).*Aux_dummy(j_LEA,:),...
                                        Z_aux(i_LEA,4).*Aux_dummy(j_LEA,:),...
                                        Z_aux(i_LEA,5).*Aux_dummy(j_LEA,:),...
                                        Z_aux(i_LEA,6).*Aux_dummy(j_LEA,:)],...
                                        beta_aux,...
                                        gamma_aux);
                                    
        COUNTERS_boots_2(i_LEA,j_LEA,i_boots) ...
                 = counterfactual_function(1,...
                                        X_aux(i_LEA,:,:),...
                                        [1,Z_aux(i_LEA,:),...
                                        Aux_dummy(i_LEA,:),...
                                        Z_aux(i_LEA,2).*Aux_dummy(i_LEA,:),...
                                        Z_aux(i_LEA,4).*Aux_dummy(i_LEA,:),...
                                        Z_aux(i_LEA,5).*Aux_dummy(i_LEA,:),...
                                        Z_aux(i_LEA,6).*Aux_dummy(i_LEA,:)],...
                                        beta_aux,...
                                        gamma_aux); 
                                    
        COUNTERS_boots(i_LEA,j_LEA,i_boots) ...                         
                 = sum(Y_aux(i_LEA,:),2).*...
                   COUNTERS_boots_1(i_LEA,j_LEA,i_boots) ./...
                   COUNTERS_boots_2(i_LEA,j_LEA,i_boots);
                                                    
        end


    end
               
end

toc;

%% 1) Quantiles for the individual effects

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

%% 2) Quantiles for total effects

lowerCI_sum ...
        = reshape(quantile(sum(COUNTERS_boots,1),alpha/2,3),...
          [1,J_LEAs]);
      
upperCI_sum ...
        = reshape(quantile(sum(COUNTERS_boots,1),1-(alpha/2),3),...
          [1,J_LEAs]);  
      
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
