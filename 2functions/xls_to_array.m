function [LEA_data, NAMES,ORI9,REGION,STATE,COUNTY,FIPS,FPLACE,AGCYTYPE,fields_names] = xls_to_array(table_LEA,variables)
%% xls_to_array: This function transforms the xls raw file into a matlab array
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [LEA_data, NAMES,ORI9,REGION,STATE,COUNTY,FIPS,fields_names] =
% xls_to_array(table_LEA,variables)
%
% If there are m variables, LEA_data is a JxTxm array
%
% INPUT:
% table_LEA: matlab created table from raw .xls file
% variables: variables to extract from the .xls
% (Field names and main identifiers are reported separately for
% manipulation)

%% 1) Create a Matlab structure with the data

ORI9s     =   unique(table_LEA.ORI9);
              %LEAs are ordered by ORI9

S         = struct([]);

for j_LEA = 1: length(ORI9s)
    
    subtable_j ...
          = table_LEA(strcmp(table_LEA.ORI9,ORI9s(j_LEA)),:);
            %subtable that selects all the variables and years available
            %for an LEA
    
    for i_variable = 1:8  
        
        S(j_LEA).(variables{i_variable}) ...
                   = subtable_j.(variables{i_variable})(1);
            %Creates one field for each of the first 7 variables in the list
            %('ORI9', 'NAME','REGION','STATENAME','COUNTYNAME','FIPS','AGCTYPE','FPLACE')     
    end
    
    for i_variable = 9:length(variables)
         
        S(j_LEA).(variables{i_variable}) ...
                   = subtable_j.(variables{i_variable})';
        
    end
    
    clear subtable_j
end  

%% 2) Create a Matlab Array for manipulation

clearvars -except S path_root


%% 2.1) Extract Field Names

fields_names = fieldnames(S);

%Extract ORI and NAMES
ORI9 = {S.ORI9}';

NAMES ...
     = {S.NAME}';
 
REGION ...
     = {S.REGION}';
 
STATE ...
     = {S.STATENAME}';
 
COUNTY ...
     = {S.COUNTYNAME}'; 
 
FIPS ...
     = {S.FIPS}';
 
FPLACE ...
     = {S.FPLACE}'; 
 
AGCYTYPE ...
     = {S.AGCYTYPE}'; 

%% 2.2) Define dimensions (T:time periods, J: Number of LEAs) 

T    = length(S(1).HOMICIDES);

J    = length(S);

fields_names ...
     = fields_names(10:end);

%% 2.3) Data array 

%Initialize
LEA_data ...
    = zeros(J,T,length(fields_names));

%Transform structure into numerical array

for i_variable = 1:length(fields_names)
   
    variable_i = extractfield(S,fields_names{i_variable});
     
    LEA_data(:,:,i_variable) ...
               = reshape(variable_i',[T,J])';         
end

%The first 2d page of LEA_data is a J x T matrix 

clear i_variable variable_i S 
%containing the time series of homicides for each LEA

%The second 2d page of LEA_data is a J x T matrix
%containing the time series of officers for each LEA 


end

