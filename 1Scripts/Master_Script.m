%% Master Script

%% 1) Call the xls file and transforms it into a matlab array
%  (approximate run time of this section: 1 minute)
%  This section calls the proprietary function:
%  xls_to_array

path_root = '../';

url       = strcat(path_root,'3Data/PD_Data_MPV_LEOKA_OKLE_CENSUS_PLUS.xlsx');

variables = {'ORI9',...
             'NAME',...             
             'REGION',...
             'STATENAME',...
             'COUNTYNAME',...
             'FIPS',...
             'AGCYTYPE',...
             'FPLACE',...
             'YEAR',...
             'HOMICIDES',...
             'OFFICERS',...
             'POPULATION',...
             'CORE_CITY',...
             'OKLE_POP',...
             'MURDERS',...
             'VIOLENT_CRIME',...
             'PROPERTY_CRIME',...
             'GUN_DEATH_RATE',...
             'BLACK_POPULATION',...
             'HISPANIC_POPULATION',...
             'BELOW_POVERTY',...
             'CENSUS_POPULATION',...
             'GARNER',...
             'LEOBR',...
             'LAND_AREA'...
%              'ASSAULTS'...
             };

only_police_departments = 1; 

addpath(strcat(path_root,'2functions'))

[LEA_data, ...
 NAMES,...
 ORI9,...
 REGION,...
 STATE,...
 COUNTY,...
 FIPS,...
 FPLACE,...
 AGCYTYPE,...
 fields_names,...
 rowsWithMissing,...
 table_missing] = create_array_file(path_root,...
                                            url,...
                                            variables,...
                                            only_police_departments);
                                        
clear variables                                        
                                        
save(strcat(path_root,'4Output/mat/LEA_data_array.mat'));

clear
%% 2) Adjust the scale of main variables

path_root = '../';

load(strcat(path_root,'4Output/mat/LEA_data_array.mat'))

OFF     = LEA_data(:,:,strcmp(fields_names,'OFFICERS'));

POP     = LEA_data(:,:,strcmp(fields_names,'POPULATION'));
          %This is a measure of the population served by each agency

OKLE_POP = LEA_data(:,:,strcmp(fields_names,'OKLE_POP'));
          %This is a measure of time-varying population 
          
MURDER ...
        = LEA_data(:,:,strcmp(fields_names,'MURDERS'));

GUN_DEATH ...
        = LEA_data(:,:,strcmp(fields_names,'GUN_DEATH_RATE'));
    
GARNER ...
       = LEA_data(:,1,strcmp(fields_names,'GARNER')); 
   
LEOBR ...
       = LEA_data(:,1,strcmp(fields_names,'LEOBR'));
   
LAND_AREA ...
       = LEA_data(:,1,strcmp(fields_names,'LAND_AREA'));  
   
% ASSAULTS ...
%        = LEA_data(:,:,strcmp(fields_names,'ASSAULTS'));
    
% Adjust the scale of Census Variables

SHARE_BLACK ...
       = 100*LEA_data(:,1,strcmp(fields_names,'BLACK_POPULATION'))...
         ./LEA_data(:,1,strcmp(fields_names,'CENSUS_POPULATION')) ;

SHARE_POVERTY ...
       = 100*LEA_data(:,1,strcmp(fields_names,'BELOW_POVERTY'))...
       ./LEA_data(:,1,strcmp(fields_names,'CENSUS_POPULATION'));    
   
%% 3) Prepare data for estimation

Y                      = LEA_data(:,:,1);

T                      = size(Y,2);

% Time-varying covariates

X(:,:,1)               = 100000*(MURDER./OKLE_POP);   %MURDER PER POPULATION PER 100,000

% X(:,:,2)               = 10*ASSAULTS./OFF;

% Time-invariant covariates

Z(:,1)                 = log(POP(:,1)./1e6);  

Z(:,2)                 = 1000*OFF(:,1)./POP(:,1);     

Z(:,3)                 = GUN_DEATH(:,1);

Z(:,4)                 = SHARE_POVERTY;

Z(:,5)                 = SHARE_BLACK;

Z(:,6)                 = GARNER; 

Z(:,7)                 = LEOBR;

Z(:,8)                 = 1e-6*LAND_AREA(:,1)./POP(:,1);

% Store variables

outcome   = {'HOMICIDES'};

time_variant_variables ...
          = {'Murder per 100k pop.'        
%           'Assaults per 10 Officers'
          };
         
time_invariant_variables ...
          = {'Log of Avg. pop. (in millions)',...             
             'Officers per 1k pop.',...
             'Gun Death Rate (\%)',...
             'Share in Poverty (\%)',...
             'Share Black (\%)',...
             'Garner',...
             'LEOBR',...
             'Land Area per 1m pop.'
              };

clearvars -except NAMES ORI9 REGION COUNTY STATE Y X Z T ...
           outcome time_variant_variables time_invariant_variables ...
           path_root

save(strcat(path_root,'4Output/mat/LEA_data_estimation_panel'));

%% 4) Dummy variable for the 10 largest cities

NAMES ...
    = cellstr(string(NAMES));

%Agencies of interest

aux ... 
    = (contains(NAMES,'PHOENIX POLICE DEPARTMENT') ) ... 
    | contains(NAMES,'LAS VEGAS METRO POLICE DEPARTMENT')...
    | (contains(NAMES,'SAN DIEGO POLICE DEPARTMENT') ) ... 
    | contains(NAMES,'LOS ANGELES POLICE DEPARTMENT')...
    | (contains(NAMES,'DALLAS POLICE DEPARTMENT'))... 
    |(contains(NAMES,'HOUSTON POLICE DEPARTMENT') )... 
    |(contains(NAMES,'PHILADELPHIA POLICE DEPARTMENT') ) ...
    | (contains(NAMES,'CHICAGO POLICE DEPT')&~contains(NAMES,'WEST'))...
    | contains(NAMES,'NEW YORK CITY POLICE DEPARTMENT')...
    | contains(NAMES,'PHILADELPHIA POLICE DEPARTMENT')...
    | contains(NAMES,'SAN ANTONIO POLICE DEPARTMENT');

aux       = aux & sum(Y,2)>1 ;
%% 5) Estimation of beta and gamma 
%  (approximate run time of this block: 2 seconds)
%  This block calls the proprietary function estimation_beta_gamma
%  10th largest cities are dropped out from estimation


addpath(strcat(path_root,'2functions'))

dummy_region ...
          = false;

[RESULTS] = estimation_beta_gamma(Y(~(aux),:),X(~(aux),:,:),Z(~(aux),:),...
                                  time_variant_variables,...
                                  time_invariant_variables,...
                                  dummy_region,...
                                  REGION(~(aux),:),...
                                  NAMES(~aux,:));

%% 6) Display results 

% Time-varying covariates: Raw Effects

K = size(X,3);

T_results_variant_raw ...
    = table(RESULTS.time_varying_covs',...
            RESULTS.betahat,...
            RESULTS.se_betahat(1:K));
        
T_results_variant_raw.Properties.VariableNames{1} ...
    = 'Time_Varying_Covariates';

T_results_variant_raw.Properties.VariableNames{2} ...
    = 'Beta';

T_results_variant_raw.Properties.VariableNames{3} ...
    = 'Standard_errors';

% display(T_results_variant_raw)

% Time-invariant covariates: Raw effects

T_results_invariant_raw ...
     = table(RESULTS.time_invariant_covs(1:end)',...
            RESULTS.gammahat(2:end),...
            RESULTS.se_gammahat(1:end)); 
        
T_results_invariant_raw.Properties.VariableNames{1} ...
    = 'Time_Invariant_Covariates';

T_results_invariant_raw.Properties.VariableNames{2} ...
    = 'Gamma';

T_results_invariant_raw.Properties.VariableNames{3} ...
    = 'Standard_errors';        

L  = size(Z,2);

display(T_results_invariant_raw)

%Produces latex for coefficients table
create_coefficient_table(RESULTS, K, L)

%% 7) Bootstrap of Beta and Gamma
% This section takes approximately 20 minutes in a MacBook Pro (Retina, 15-inch, Early 2013)
% (running on MacOs Catalina)


nboot = 1000;

alpha = .1;

[RESULTS_boot] = efron_bootstrap_se(nboot, alpha, Y(~aux,:),X(~aux,:,:),Z(~aux,:),time_variant_variables,time_invariant_variables);                              

bootstrap_create_coefficient_table(RESULTS_boot, K, L)
%% 8) Counterfactuals for Top 10 largest Agencies: Unobservables
%Confidence Intervals are generated by sampling from the asymptotic
%distribution of beta and gamma

NAMES_aux = NAMES(aux);
 
REGION_aux = REGION(aux);

X_aux     = X(aux,:,:);
 
Y_aux     = Y(aux,:);

if dummy_region == false 
 
Z_aux ...
          = [ones(sum(aux,1),1),Z(aux,:)]; %includes a constant
      
else
   
Z_aux = create_dummy_region(REGION_aux,NAMES_aux,Z(aux,:));
  
Z_aux = [ones(sum(aux,1),1),Z_aux]; %includes a constant and dummies
 
end

confidence_level = .90;
      
[NAMES_aux_sorted, ~, CIs,CIs_sum,NAIVE_COUNTER,~] = master_counter(NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,confidence_level);  
                                          
display(NAMES_aux_sorted)

display(CIs)

display(CIs_sum) 

% Produces latex for counterfactuals table
create_table_unobs(NAMES_aux_sorted,CIs,CIs_sum)

%% 9) Bootstrapped Counterfactuals for Top 10 largest Agencies: Unobservables
%Confidence intervals are generated by using the bootstrap (resampling agencies)

[NAMES_aux_sorted, ~, CIs,CIs_sum, ~] = bootstrap_master_counter(nboot,NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,RESULTS_boot,confidence_level);

bootstrap_create_table_unobs(NAMES_aux_sorted,CIs,CIs_sum)

%% 10) Counterfactuals for Top 10 largest Agencies: Observables
%Confidence Intervals are generated by sampling from the asymptotic
%distribution of beta and gamma

policy_vars = [3 4 5 6 7];
%Officersper1K, gundeathrate, share of poverty, Garner, LEOBR

[NAMES_aux_sorted,...
~,...
CIs,...
CIs_sum,...
NAIVE_counter,...
~] = master_counter_observables(NAMES_aux,...
                                REGION_aux,...
                                Y_aux,X_aux,Z_aux,...
                                policy_vars,...
                                RESULTS,...
                                confidence_level);  
                                          
display(NAMES_aux_sorted)
 
display(CIs)
  
display(CIs_sum) 

%Produces latex for counterfactuals table
create_table_obs(NAMES_aux_sorted, CIs, CIs_sum)

%% 11) Bootstrapped Counterfactuals for Top 10 largest Agencies: Observables
%Confidence intervals are generated by using the bootstrap (resampling agencies)

[NAMES_aux_sorted, ~, CIs,CIs_sum, post_mean_hat_sorted] = bootstrap_master_counter_observables(policy_vars,nboot,NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,RESULTS_boot,confidence_level);

display(NAMES_aux_sorted)

display(CIs)
  
display(CIs_sum) 

bootstrap_create_table_obs(NAMES_aux_sorted,CIs,CIs_sum)

%% 12) Counterfactuals for Top 10 largest Agencies: Observables and Unobservables
%Confidence Intervals are generated by sampling from the asymptotic
%distribution of beta and gamma

policy_vars = [3 4 5 6 7];

[NAMES_aux_sorted, ~,CIs,CIs_sum,~] = master_counter_unobs_obs(NAMES_aux,...
                                REGION_aux,...
                                Y_aux,X_aux,Z_aux,...
                                policy_vars,...
                                RESULTS,...
                                confidence_level);  

display(NAMES_aux_sorted)

display(CIs)
  
display(CIs_sum)                             
                            
%Produces latex for counterfactuals table
create_table_unobs_obs(NAMES_aux_sorted, CIs, CIs_sum)
%% 13) Bootstrapped Counterfactuals for Top 10 largest Agencies: Observables and Unobservables
%Confidence intervals are generated by using the bootstrap (resampling agencies)

policy_vars = [3 4 5 6 7];

[NAMES_aux_sorted, REGION_aux_sorted,CIs,CIs_sum,alpha_hat_sorted] = bootstrap_master_counter_unobs_obs(nboot,NAMES_aux,...
                                REGION_aux,...
                                Y_aux,X_aux,Z_aux,...
                                policy_vars,...
                                RESULTS,...
                                RESULTS_boot,...
                                confidence_level);  

display(NAMES_aux_sorted)

display(CIs)
  
display(CIs_sum) 
                            
% Produces latex for counterfactuals table
bootstrap_create_table_unobs_obs(NAMES_aux_sorted, CIs, CIs_sum)


%% 16) Assets
%Generates different values that appear in the paper (number of agencies,
%number of homicides, etc).

create_assets(Y,CIs,CIs_sum,NAIVE_COUNTER,RESULTS)
