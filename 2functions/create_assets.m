function create_assets(Y,CIs_unobs,CIs_sum_unobs,CIs_obs,CIs_sum_obs,CIs_unob_obs,CIs_sum_unobs_obs,RESULTS)
%% CREATE_ASSETS: This function creates and saves Latex assets 
%
% Author: José Luis Montiel Olea. Last Revised: June 20th, 2021
%
% Built using MATLAB Version: 9.7.0.1296695 (R2019b) Update 4
%
% Syntax: create_assets(Y,CIs_unobs,CIs_sum_unobs,CIs_obs,CIs_sum_obs,CIs_unob_obs,CIs_sum_unobs_obs,RESULTS) 
%
% INPUT:
% Y: matrix of size J x T
% CIs: Confidence intervals for counterfactual values
% CIs_sum: Confidence intervals summed across departments for a given a_j
% RESULTS: Result of beta_gamma_estimation
% 
%% 

dfile = '../4Output/tex/Assets.tex';
if exist(dfile,'file'); delete(dfile); end 
diary(dfile)
diary on 

%finding the number of agencies with positive homicides
array_homicides = sum(Y');
r = 0;
for i = 1:length(array_homicides)
    if array_homicides(i) > 0
        r = r + 1;
    end
end 

%extracting homicides per city from CIs
ph = split(CIs_obs(1,1),','); ph2 = ph(2,1); ph3 = split(ph2,']'); ph4 = ph3(1,1);
lv = split(CIs_obs(2,2),','); lv2 = lv(2,1); lv3 = split(lv2,']'); lv4 = lv3(1,1);
dl = split(CIs_obs(3,3),','); dl2 = dl(2,1); dl3 = split(dl2,']'); dl4 = dl3(1,1);
la = split(CIs_obs(4,4),','); la2 = la(2,1); la3 = split(la2,']'); la4 = la3(1,1);
sa = split(CIs_obs(5,5),','); sa2 = sa(2,1); sa3 = split(sa2,']'); sa4 = sa3(1,1);
hu = split(CIs_obs(6,6),','); hu2 = hu(2,1); hu3 = split(hu2,']'); hu4 = hu3(1,1);
sd = split(CIs_obs(7,7),','); sd2 = sd(2,1); sd3 = split(sd2,']'); sd4 = sd3(1,1);
ch = split(CIs_obs(8,8),','); ch2 = ch(2,1); ch3 = split(ch2,']'); ch4 = ch3(1,1);
pl = split(CIs_obs(9,9),','); pl2 = pl(2,1); pl3 = split(pl2,']'); pl4 = pl3(1,1);
ny = split(CIs_obs(10,10),','); ny2 = ny(2,1); ny3 = split(ny2,']'); ny4 = ny3(1,1);

tt = str2double(ph4{1,1}) + str2double(lv4{1,1}) + str2double(dl4{1,1}) + str2double(la4{1,1}) + str2double(sa4{1,1}) + str2double(hu4{1,1}) + str2double(sd4{1,1}) + str2double(ch4{1,1}) + str2double(pl4{1,1}) + str2double(ny4{1,1});
nyt = split(CIs_sum_unobs_obs(1,10),','); nyt2 = nyt(2,1); nyt3 = split(nyt2,']'); nyt4 = nyt3(1,1);

%defining the assets
totalagencies = string(size(Y,1));
totalincidentsten = string(tt);
peopleserved = '1.35 million';
totalhomicides = '3820';
percentagehomicides = string(round(((tt/str2double(totalhomicides))*100)));
newyorktotalunobs = string(CIs_sum_unobs(1,10));
phoenixtotalunobs = string(CIs_sum_unobs(1,1));
phoenix = string(ph4);
phoenixnewyorkunobs = string(CIs_unobs(1,10));
newyork = string(ny4);
newyorkphoenixunobs = string(CIs_unobs(10,1));
officerscoeff = string(round((RESULTS.gammahat(3)*100),2));
gundeathratecoeff = string(round((RESULTS.gammahat(4)*100),1));
shareinpoverty = string(round((RESULTS.gammahat(5)*100)));
newyorktotalunobsobs = string(CIs_sum_unobs_obs(1,10));
totalvictims = '6574';
ondutyhomicides = '5885';
policedepartments = '13397';
leokadepartments = '12012';
leokahomicides = '3799';
okledepartments = '7720';
oklehomicides = '3484';
povertydepartments = '7667';
povertyhomicies = '3477';
positivehomicides = string(r);
finalhomicides = '3476';
officerscoeffrounded = string(round(RESULTS.gammahat(3),3));
povertypercentage = '4.68';
avgyearlyencounters = string(round(str2double(finalhomicides)/6));
reductionpoverty = string(round(str2double(avgyearlyencounters)*(str2double(povertypercentage)/100)));
newyorkreduction = string(str2double(totalincidentsten) - str2double(nyt4));
philadelphia = string(pl4);
philanewyorkunobs = string(CIs_unobs(4,10));
philanewyorkobs = string(CIs_obs(9,10));
newyorktotalobs = string(CIs_sum_obs(1,10));

array1 = [totalagencies,...
totalincidentsten,...
peopleserved,...
totalhomicides,...
percentagehomicides,...
newyorktotalunobs,...
phoenixtotalunobs,...
phoenix,...
phoenixnewyorkunobs,...
newyork,...
newyorkphoenixunobs,...
officerscoeff,...
gundeathratecoeff,...
shareinpoverty,...
newyorktotalunobsobs,...
totalvictims,...
ondutyhomicides,...
policedepartments,...
leokadepartments,...
leokahomicides,...
okledepartments,...
oklehomicides,...
povertydepartments,...
povertyhomicies,...
positivehomicides,...
finalhomicides,...
officerscoeffrounded,...
povertypercentage,...
avgyearlyencounters,...
reductionpoverty,...
newyorkreduction,...
philadelphia,...
philanewyorkunobs,...
philanewyorkobs,...
newyorktotalobs];
     
array2 = ["totalagencies",...
         "totalincidentsten",...
         "peopleserved",...
         "totalhomicides",...
         "percentagehomicides",...
         "newyorktotalunobs",...
         "phoenixtotalunobs",...
         "phoenix",...
         "phoenixnewyorkunobs",...
         "newyork",...
         "newyorkphoenixunobs",...
         "officerscoeff",...
         "gundeathratecoeff",...
         "shareinpoverty",...
         "newyorktotalunobsobs",...
         "totalvictims",...
         "ondutyhomicides",...
         "policedepartments",...
         "leokadepartments",...
         "leokahomicides",...
         "okledepartments",...
         "oklehomicides",...
         "povertydepartments",...
         "povertyhomicies",...
         "positivehomicides",...
         "finalhomicides",...
         "officerscoeffrounded",...
         "povertypercentage",...
         "avgyearlyencounters",...
         "reductionpoverty",...
         "newyorkreduction",...
         "philadelphia",...
         "philanewyorkunobs",...
         "philanewyorkobs",...
         "newyorktotalobs"];

%saving the assets to a .tex file     
for i = 1:length(array1)
    path = strcat('../4Output/tex/Assets.tex');
    fid = fopen(path,'at');
    str = strcat("\\newcommand{\\",array2(i),"}{",array1(i),"}");
    fprintf(fid, str);
    fprintf(fid,'\n');
    fclose(fid);
end
    

     
    
     
        


