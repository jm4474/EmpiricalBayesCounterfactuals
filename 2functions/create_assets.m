function create_assets(Y,CIs,CIs_sum,NAIVE_counter)
%% CREATE_ASSETS: This function creates and saves Latex assets 
%
% Syntax: create_assets(Y,CIs,CIs_sum,NAIVE_counter) 
%
% INPUT:
% Y:   
% CIs: 
% CIs_sum:
% NAIVE_counter:
% 
%% 

totalagencies = string(size(Y,1));
totaldeadlyforceten = string(trace(NAIVE_counter));
phoenix = string(NAIVE_counter(1,1));
newyork = string(NAIVE_counter(10,10));
phoenixtotalcounterfac = CIs_sum(1,1);
newyorktotalcounterfac = CIs_sum(1,10);
phoenixnewyorkcounterfac = CIs(1,10);
newyorkphoenixcounterfac = CIs(10,1);

array1 = [totalagencies,...
         totaldeadlyforceten,...
         phoenix,...
         newyork,...
         phoenixtotalcounterfac,...
         newyorktotalcounterfac,...
         phoenixnewyorkcounterfac,...
         newyorkphoenixcounterfac];
     
array2 = ["totalagencies",...
         "totaldeadlyforceten",...
         "phoenix",...
         "newyork",...
         "phoenixtotalcounterfac",...
         "newyorktotalcounterfac",...
         "phoenixnewyorkcounterfac",...
         "newyorkphoenixcounterfac"];
     
for i = 1:length(array1)
    path = strcat('../4Output/assets/',array2(i),'.tex');
    fid = fopen(path,'wt');
    str = strcat("\\newcommand{\\",array2(i),"}{",array1(i),"}");
    fprintf(fid, str);
    fclose(fid);
end 
    

     
    
     
        


