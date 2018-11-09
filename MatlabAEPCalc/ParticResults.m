% Nick Baker
% Created for IEA Task 37 Wind Farm Layout Optimization Case Study 1
% To assess participant results
% Created 30 Oct 18

clear, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

fname_turb_loc = 'iea37-test9.yaml';
farmSize = 0;   % 0 = 16 turbines, 1 = 36 turbines, 2 = 64 turbines
switch(farmSize)
    case 0
        farm_rad = 900;
        plot_dimen = 1200;
        nNumRtrs = 9;
    case 1
        farm_rad = 1300;
        plot_dimen = 1600;
        nNumRtrs = 16;
    case 2
        farm_rad = 2000;
        plot_dimen = 2500;
        nNumRtrs = 36;
    case 3
        farm_rad = 3000;
        plot_dimen = 3500;
        nNumRtrs = 64;
    otherwise
        error('Variable "FarmSize" not initilized properly');
end

% Get turbine data from .yaml
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML(fname_turb_loc);
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);

% Plot their farm
color_num = 2;  % 0 = blue, 1 = red, 2 = yellow
plotFarm(turb_coords, turb_diam, farm_rad, plot_dimen, color_num)