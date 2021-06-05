function [table] = bootstrap_create_coefficient_table(RESULTS_boot, K, L)
%% BOOTSTRAP_CREATE_COEFFICIENT_TABLE: This function creates Latex that can be used to generate a table of coefficients beta and gamma
%
% Author: José Luis Montiel Olea. Last Revised: May 29th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [table] = bootstrap_create_coefficient_table(RESULTS_boot, K, L)
%
% INPUT:
% RESULTS: Result of beta_gamma_estimation
% K: size(X,3)
% L: size(Z,2)
%% 
diary '../4Output/tex/Bootstrap_Table_Coefficients.txt'
fprintf('\\begin{table}[ht]\n\\centering\n\\caption{Coefficient Estimates} \\label{table:coefficients}\n\\begin{tabular}{lcc}\n\\hline \\hline\n& Confidence Intervals & Standard Errors \\\\ \n\\hline\n')
for i = 1:K
    T = strcat(RESULTS_boot.time_variant_variables(i), ' & ',strcat('[',string(round(RESULTS_boot.betas_CI(i,1),2)),',',string(round(RESULTS_boot.betas_CI(i,2),2)),']') ,'&', string(round(RESULTS_boot.betas_SE(i), 3)), '\\');
    fprintf('%s\n', T)
end 
for i = 1:L;
    T = strcat(RESULTS_boot.time_invariant_covs(i), ' & ',strcat('[',string(round(RESULTS_boot.gammas_CI(i+1,1),2)),',',string(round(RESULTS_boot.gammas_CI(i+1,2),2)),']'),'&', string(round(RESULTS_boot.gammas_SE(i+1), 3)), '\\');
    fprintf('%s\n', T)
end 
fprintf('\\hline \\end{tabular} \\end{table}')
diary off
end
