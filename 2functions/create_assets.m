function create_assets(Y,CIs,CIs_sum,NAIVE_counter,RESULTS)
%% CREATE_ASSETS: This function creates and saves Latex assets 
%
% Author: José Luis Montiel Olea. Last Revised: June 15th, 2021
%
% Built using MATLAB Version: 9.7.0.1296695 (R2019b) Update 4
%
% Syntax: create_assets(Y,CIs,CIs_sum,NAIVE_counter) 
%
% INPUT:
% Y: matrix of size J x T
% CIs: Confidence intervals for counterfactual values
% CIs_sum: Confidence intervals summed across departments for a given a_j
% NAIVE_counter:
%RESULTS:
% 
%% 

total_agencies = string(size(Y,1));
total_deadly_force_ten = string(trace(NAIVE_counter));
phoenix = string(NAIVE_counter(1,1));
newyork = string(NAIVE_counter(10,10));
phoenix_total_counterfac = CIs_sum(1,1);
newyork_total_counterfac = CIs_sum(1,10);
phoenix_newyork_counterfac = CIs(1,10);
newyork_phoenix_counterfac = CIs(10,1);
gun_death_rate_coeff = string((RESULTS.gammahat(4)*100));
officers_per_k_coeff = string((RESULTS.gammahat(3)*100));


array1 = [total_agencies,...
         total_deadly_force_ten,...
         phoenix,...
         newyork,...
         phoenix_total_counterfac,...
         newyork_total_counterfac,...
         phoenix_newyork_counterfac,...
         newyork_phoenix_counterfac,...
         gun_death_rate_coeff,...
         officers_per_k_coeff];
     
array2 = ["total_agencies",...
         "total_deadly_force_ten",...
         "phoenix",...
         "newyork",...
         "phoenix_total_counterfac",...
         "newyork_total_counterfac",...
         "phoenix_newyork_counterfac",...
         "newyork_phoenix_counterfac",...
         "gun_death_rate_coeff",...
         "officers_per_k_coeff"];
     
for i = 1:length(array1)
    path = strcat('../4Output/assets/',array2(i),'.tex');
    fid = fopen(path,'wt');
    str = strcat("\\newcommand{\\",array2(i),"}{",array1(i),"}");
    fprintf(fid, str);
    fclose(fid);
end 
    

     
    
     
        


