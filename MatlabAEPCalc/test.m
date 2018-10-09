
clear all, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

% Get turbine data from .yaml
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex16.yaml');
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);
nNumRtrs = length(turb_coords.x);   % Pulls the number of turbines by how many x-coordinates we have.
[nNumPairs,~] = size(combnk(1:nNumRtrs,2));   % Calculates number of unique pairs for spacing constraints using combinatorics 

nConstIndex = 0;            % Initialize our index array
c1 = zeros(nNumRtrs,1);     % Initialize constraint array, for how many unique pairs we have, plus one spot for each turbine
for i = 1:nNumRtrs
    % Farm boundary constraint
    c1(i) = -(turb_coords.x(i)^2) - (turb_coords.y(i)^2) + (turb_diam/2)^2;
end



(c1 == c2)
