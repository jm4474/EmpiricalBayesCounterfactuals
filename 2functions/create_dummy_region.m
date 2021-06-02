function [Z] = create_dummy_region(REGION,NAMES,Z)
%% create_dummy_region: This function creates dummy variables for each region 
%
% Author: José Luis Montiel Olea. Last Revised: May 24th, 2021
%
% Built using MATLAB Version: 9.10.0.1538726 (R2021a) Prerelease
%
% Syntax: [Z] = create_dummy_region(REGION,NAMES,Z)
%
% INPUT:
% REGION: Array of regions for the data
% NAMES: Array of names for the data
% Z: Array of size  J x L

%% 
REGION ...
    = cellstr(string(REGION));

NAMES ...
    = cellstr(string(NAMES));

REGION_UNIQUE ...
    = unique(REGION);

Number_regions ...
    = size(REGION_UNIQUE,1);

dummy_region ...
    = zeros(size(NAMES,1),Number_regions-1);

for i_region = 1:Number_regions-1
   
    dummy_region(:,i_region) ...
             = strcmp(REGION_UNIQUE(i_region),REGION);
         
end    

Z = [Z,...     
     dummy_region(:,1),...     
     dummy_region(:,2),...
     dummy_region(:,3),...
     dummy_region.*Z(:,2),...
     dummy_region.*Z(:,4),...
     dummy_region.*Z(:,5),...
     dummy_region.*Z(:,6)];


end

