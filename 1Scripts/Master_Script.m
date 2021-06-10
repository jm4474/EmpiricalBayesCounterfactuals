%% Master Script

%% 1) Adjust the scale of main variables

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
    
% Adjust the scale of Census Variables

SHARE_BLACK ...
       = 100*LEA_data(:,1,strcmp(fields_names,'BLACK_POPULATION'))...
         ./LEA_data(:,1,strcmp(fields_names,'CENSUS_POPULATION')) ;

SHARE_POVERTY ...
       = 100*LEA_data(:,1,strcmp(fields_names,'BELOW_POVERTY'))...
       ./LEA_data(:,1,strcmp(fields_names,'CENSUS_POPULATION'));    
   
%% 2) Prepare data for estimation

Y                      = LEA_data(:,:,1);

T                      = size(Y,2);

% Time-varying covariates

X(:,:,1)               = 100000*(MURDER./OKLE_POP);   %MURDER PER POPULATION PER 100,000

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
          = {'Murder per 100k pop.'};
         
time_invariant_variables ...
          = {'Log of Avg. pop. (in millions)',...             
             'Officers per 1k pop.',...
             'Gun Death Rate (\%)',...
             'Share in Poverty (\%)',...
             'Share Black (\%)',...
             'Garner',...
             'LEOBR',...
             'Land Area per pop.'
              };

clearvars -except NAMES ORI9 REGION COUNTY STATE Y X Z T ...
           outcome time_variant_variables time_invariant_variables ...
           path_root

save(strcat(path_root,'4Output/mat/LEA_data_estimation_panel'));

%% 3) Dummy variable for the 10 largest cities

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

%% 4) Estimation of beta and gamma 
%  (approximate run time of this block: 2 seconds)
%  This block calls the proprietary function estimation_beta_gamma
%  10th largest cities are dropped out from estimation

addpath(strcat(path_root,'2functions'))

dummy_region ...
          = false;
      
[RESULTS] = estimation_beta_gamma(Y(~aux,:),X(~aux,:,:),Z(~aux,:),...
                                  time_variant_variables,...
                                  time_invariant_variables,...
                                  dummy_region,...
                                  REGION(~aux,:),...
                                  NAMES(~aux,:));

%% 5) Display results 

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

% display(T_results_invariant_raw)

%Produces latex for coefficients table
create_coefficient_table(RESULTS, K, L)

%% 6) Display results: Percentages

K                      = size(X,3);

T_results_variant ...
    = table(RESULTS.time_varying_covs',...
            100*(exp(RESULTS.betahat)-1),...
            100*abs(exp(RESULTS.betahat)).*RESULTS.se_betahat(1:K));
        
T_results_variant.Properties.VariableNames{1} ...
    = 'Time_Varying_Covariates';

T_results_variant.Properties.VariableNames{2} ...
    = 'Estimated_effect';

T_results_variant.Properties.VariableNames{3} ...
    = 'Standard_errors';

display(T_results_variant)

% Percent Effect of changing population 100x%

 x= .1; 

T_results_invariant ...
     = table(RESULTS.time_invariant_covs(1)',...
             100*((1+x).^(RESULTS.gammahat(2))-1),...
             abs(100*log(1+x).*exp(log(1+x).*RESULTS.gammahat(2))) ...
             .*RESULTS.se_gammahat(1)); 
        
% Percent Effect of changing covariates in one unit
        
T_results_invariant ....
    = [T_results_invariant; table(RESULTS.time_invariant_covs(2:end)',...
            100*(exp(RESULTS.gammahat(3:end))-1),...
            100*exp(RESULTS.gammahat(3:end)) ...
            .*RESULTS.se_gammahat(2:end))];   
        
        
T_results_invariant.Properties.VariableNames{1} ...
    = 'Time_Invariant_Covariates';

T_results_invariant.Properties.VariableNames{2} ...
    = 'Estimated_effect';

%(10% for pop and crime, and 1 more off per pop)

T_results_invariant.Properties.VariableNames{3} ...
    = 'Standard_errors';  

display(T_results_invariant)

%% 7) Bootstrap of Beta and Gamma

nboot = 100;
alpha = .1;

[RESULTS_boot] = efron_bootstrap_se(nboot, alpha, Y(~aux,:),X(~aux,:,:),Z(~aux,:),time_variant_variables,time_invariant_variables);                              

bootstrap_create_coefficient_table(RESULTS_boot, K, L)
%% 8) Counterfactuals for Top 10 largest Agencies: Unobservables

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
      
[NAMES_aux_sorted, REGION_aux_sorted, CIs,CIs_sum,NAIVE_COUNTER,post_mean_hat_sorted] = master_counter(NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,confidence_level);  
                                          
% display(NAMES_aux_sorted)
%
% display(CIs)
%
% display(CIs_sum) 

% Produces latex for counterfactuals table
create_table_unobs(NAMES_aux_sorted,CIs,CIs_sum)

%% 9) Bootstrapped Counterfactuals for Top 10 largest Agencies: Unobservables

[NAMES_aux_sorted, REGION_aux_sorted, CIs,CIs_sum, post_mean_hat_sorted] = bootstrap_master_counter(nboot,NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,RESULTS_boot,confidence_level);

bootstrap_create_table_unobs(NAMES_aux_sorted,CIs,CIs_sum)

%% 10) Counterfactuals for Top 10 largest Agencies: Observables

policy_vars = [3 4 5 6 7];

[NAMES_aux_sorted,...
REGION_aux_sorted,...
CIs,...
CIs_sum,...
NAIVE_counter,...
alpha_hat_sorted] = master_counter_observables(NAMES_aux,...
                                REGION_aux,...
                                Y_aux,X_aux,Z_aux,...
                                policy_vars,...
                                RESULTS,...
                                confidence_level);  
                                          
% display(NAMES_aux_sorted)
% 
% display(CIs)
%  
% display(CIs_sum) 

% Produces latex for counterfactuals table
create_table_obs(NAMES_aux_sorted, CIs, CIs_sum)

%% 11) Bootstrapped Counterfactuals for Top 10 largest Agencies: Observables

[NAMES_aux_sorted, REGION_aux_sorted, CIs,CIs_sum, post_mean_hat_sorted] = bootstrap_master_counter_observables(policy_vars,nboot,NAMES_aux,REGION_aux,Y_aux,X_aux,Z_aux,RESULTS,RESULTS_boot,confidence_level);

bootstrap_create_table_obs(NAMES_aux_sorted,CIs,CIs_sum)

%% 12) Counterfactuals by regions

%This section requires the creation of artificial regional LEAs

NAMES_synthetic = unique(cellstr(string(REGION)));

load('../LEA_artificial_regional_agencies.mat')

if dummy_region == false 
    
Z_region ...
             = [ones(size(Z_region,1),1),Z_region]; %includes a constant

else

Z_region     = create_dummy_region({'EAST';'MIDWEST';'SOUTH';'WEST'},...
                                   {'EAST';'MIDWEST';'SOUTH';'WEST'},...
                                   Z_region);  
Z_region ...
          = [ones(size(Z_region,1),1),Z_region]; %includes a constant and dummies
end

[NAMES_synthetic_aux_sorted,...
REGION_synthetic_aux_sorted,...
CIs_synthetic,...
CIs_synthetic_sum,...
NAIVE_counter_synthtetic,...
alpha_hat_synthetic_sorted] = master_counter(NAMES_synthetic,...
                                NAMES_synthetic,...
                                Y_region,X_region,Z_region,...
                                RESULTS,...
                                confidence_level);  
display(NAMES_synthetic_aux_sorted)

display(CIs_synthetic)

display(CIs_synthetic_sum)  

%% 13) Assets

create_assets(Y,CIs,CIs_sum,NAIVE_counter)
