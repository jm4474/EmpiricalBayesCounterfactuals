function [table] = create_coefficient_table(RESULTS, K, L)
%% CREATE_COEFFICIENT_TABLE: This function creates Latex that can be used to generate a table of coefficients beta and gamma
%
% Author: José Luis Montiel Olea. Last Revised: May 29th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [table] = create_coefficient_table(RESULTS, K, L)
%
% INPUT:
% RESULTS: Result of beta_gamma_estimation
% K: size(X,3)
% L: size(Z,2)
%% 
diary '../4Output/tex/Table_Coefficients.txt'
fprintf('\\begin{table}[ht]\n\\centering\n\\caption{Coefficient Estimates} \\label{table:coefficients}\n\\begin{tabular}{lcc}\n\\hline \\hline\n& Coefficient & Standard Errors \\\\ \n\\hline\n')
for i = 1:K
    T = strcat(RESULTS.time_varying_covs(i), ' & ',string(round(RESULTS.betahat(i), 3)),'&', string(round(RESULTS.se_betahat(i), 3)), '\\');
    fprintf('%s\n', T)
end 
for i = 1:L;
    T = strcat(RESULTS.time_invariant_covs(i), ' & ',string(round(RESULTS.gammahat(i+1), 3)),'&', string(round(RESULTS.se_gammahat(i), 3)), '\\');
    fprintf('%s\n', T)
end 
fprintf('\\hline \\end{tabular} \\end{table}')
diary off
end
