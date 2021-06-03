function [table] = create_table_unobs(NAMES_aux_sorted,CIs,CIs_sum)
%% CREATE_TABLE_UNOBS: This function creates Latex that can be used to generate a table of counterfactual values of lethal encounters
%
% Author: José Luis Montiel Olea. Last Revised: May 29th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [table] = create_table_unobs(NAMES_aux_sorted,CIs,CIs_sum)
%
% INPUT:
% NAMES_aux_sorted: Names of police departments sorted by parameter a_j
% CIs: Confidence intervals for counterfactual values
% CIs_Sum: Confidence intervals summed across departments for a given a_j
%%
diary '../4Output/tex/Table_Counterfactual_Unobs.txt'

fprintf('\\begin{table}[h!]\\centering\\scriptsize\\tabcolsep=0.1cm\\caption{Counterfactual police homicides for 2013-2018: Unobserved Covariates} \\label{table:counterfactual_unobs}\\begin{tabular}{lrrrrrrrrrr}\\hline \\hline')
T = '';
for i = 1:size(NAMES_aux_sorted,1)
    str = string(NAMES_aux_sorted(i,1));
    str = lower(str);
    expression = '(^|[\. ])\s*.';
    replace = '${upper($0)}';
    str = regexprep(str,expression,replace);
    str = strsplit(str);
    if str(1) == 'Los' | str(1) == 'Las' | str(1) == 'San'| str(1) == 'New'
        T = append(T,strcat('&',strcat(str(1),{' '},str(2))));
    else 
        T = append(T,strcat('&',str(1)));
    end
end
T = string(strcat(T, '\\'));
fprintf('%s\n', T)
E = '\hline';
fprintf('%s\n', E)
for i = 1:size(NAMES_aux_sorted,1)
    T = '';
    str = string(NAMES_aux_sorted(i,1));
    str = lower(str);
    expression = '(^|[\. ])\s*.';
    replace = '${upper($0)}';
    str = regexprep(str,expression,replace);
    str = strsplit(str);
    if str(1) == 'Los' | str(1) == 'Las' | str(1) == 'San'| str(1) == 'New'
        T = append(T,strcat(str(1),{' '},str(2)));
    else 
        T = append(T, str(1));
    end
    for j = 1:size(NAMES_aux_sorted,1)
        if i == j
            T = append(T,strcat('&','\textbf{',strip(strip(eraseBetween(strip(string(CIs(i,j)),'['), ',', ']'),']'),','),'}'));
        else 
            T = append(T, strcat('&',string(CIs(i,j))));
        end 
    end
    T = strcat(T, '\\');
    fprintf('%s\n', T)
end
fprintf('%s\n', E)
T = 'Totals';
for i = 1:size(CIs_sum,2)
    T = append(T, strcat('&',string(CIs_sum(1,i))));
end
T = strcat(T, '\\');
fprintf('%s\n', T)
fprintf('\\hline\\end{tabular}\\begin{center}\\begin{minipage}{1.05\\textwidth} %% choose width suitably \n {\\footnotesize {\\schape Note}: Diagonal entries are observed lethal encounters. Off-diagonal entries are 90\\%% confidence intervals for counterfactual values of lethal encounters, obtained by replacing the posterior expectation of $\\alpha_j$ for the agency in the row with that of the agency in the column, while assuming no change in observed covariates. Agencies are listed in decreasing order by their estimated value of $\\alpha_j$.} \n \\end{minipage} \n \\end{center} \n \\end{table} \n')

diary off 
end

