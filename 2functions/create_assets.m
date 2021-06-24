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
% RESULTS: Result of beta_gamma_estimation
% NAIVE_counter:
% 
%% 

array_homicides = sum(Y');
r = 0;
for i = 1:length(array_homicides)
    if array_homicides(i) > 0
        r = r + 1;
    end
end 

positivehomicides = string(r);
totalagencies = string(size(Y,1));
totalincidentsten = string(trace(NAIVE_counter));
phoenix = string(NAIVE_counter(1,1));
newyork = string(NAIVE_counter(10,10));
phoenixtotalcounterfac = CIs_sum(1,1);
newyorktotalcounterfac = CIs_sum(1,10);
phoenixnewyorkcounterfac = CIs(1,10);
newyorkphoenixcounterfac = CIs(10,1);
gundeathratecoeff = string((RESULTS.gammahat(4)*100));
officerscoeff = string((RESULTS.gammahat(3)*100));


array1 = [positivehomicides,...
         totalagencies,...
         totalincidentsten,...
         phoenix,...
         newyork,...
         phoenixtotalcounterfac,...
         newyorktotalcounterfac,...
         phoenixnewyorkcounterfac,...
         newyorkphoenixcounterfac,...
         gundeathratecoeff,...
         officerscoeff];
     
array2 = ["positivehomicides",...
         "totalagencies",...
         "totalincidentsten",...
         "phoenix",...
         "newyork",...
         "phoenixtotalcounterfac",...
         "newyorktotalcounterfac",...
         "phoenixnewyorkcounterfac",...
         "newyorkphoenixcounterfac",...
         "gundeathratecoeff",...
         "officerscoeff"];
     
for i = 1:length(array1)
    path = strcat(''../4Output/tex/Assets.tex');
    fid = fopen(path,'at');
    str = strcat("\\newcommand{\\",array2(i),"}{",array1(i),"}");
    fprintf(fid, str);
    fprintf(fid,'\n');
    fclose(fid);
end 
    

     
    
     
        


