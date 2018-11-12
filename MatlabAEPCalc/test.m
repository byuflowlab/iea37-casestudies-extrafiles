
clear, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));
addpath(genpath('/Users/nbaker/GitHub/iea37-casestudies-extrafiles/MatlabAEPCalc'));

farmSize = 1;   % 0 = 9 turbines, 1 = 16 turbs, 2 = 36 turbs, 3 = 64 turbs
switch(farmSize)
    case 0
        fname_turb_loc = 'iea37-par1-cc1.yaml';
        fig_name = 'iea37-turbs-prakash16.pdf';
        farm_rad = 900;
        plot_dimen = 1200;
    case 1
        fname_turb_loc = 'iea37-prakash-opt16.yaml';
        fig_name = 'iea37-turbs-prakash16.pdf';
        farm_rad = 1300;
        plot_dimen = 1600;
    case 2
        fname_turb_loc = 'iea37-prakash-opt36.yaml';
        fig_name = 'iea37-turbs-prakash36.pdf';
        farm_rad = 2000;
        plot_dimen = 2500;
    case 3
        fname_turb_loc = 'iea37-prakash-opt64.yaml';
        fig_name = 'iea37-turbs-prakash64.pdf';
        farm_rad = 3000;
        plot_dimen = 3500;
    otherwise
        error('Variable "FarmSize" not initilized properly');
end

% Get turbine data from .yaml
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML(fname_turb_loc);
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);
nNumRtrs = length(turb_coords.x);   % Pulls the number of turbines by how many x-coordinates we have.

% Get AEP data from turb locations
binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
AEP = sum(binned_AEP)

color_num = 6;  % 0 = blue, 1 = red, 2 = yellow, 3 = purple, 4 = green
plotFarm(turb_coords, turb_diam, farm_rad, plot_dimen, color_num)
saveas(gcf,fig_name)
