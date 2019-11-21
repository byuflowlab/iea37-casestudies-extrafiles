% Nick Baker
% iea37 cs3&4
% test file display regions due to concavities
clear, close all;
addpath('../startup-files/')
addpath('../figure-code/')
%addpath('turb-layout/')
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

%-- Read in YAML files --%
%[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex-opt4.yaml');
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex-opt3.yaml');
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_dir_freq, wind_speed, wind_speed_freq] = getWindroseYaml(fname_wr);
num_speed_bins = length(wind_speed);
%[IIIa, IIIb, IVa, IVb, IVc] = getBorsBoundariesYaml('iea37-boundary-cs4.yaml');
[IIIa] = getBorsBoundaryYaml('iea37-boundary-cs3.yaml');

%-- Show Farm --%
hold on
%plotBorsseleBoundaries(IIIa, IIIb, IVa, IVb, IVc)
plotBorsseleBoundary(IIIa)
%plotBorsseleSlices(IIIa)
indexPts = [1,7,9,10];
%plotBorsseleMapping(IIIa, indexPts)
plotBorsseleSlices(IIIa)
%plotBorselleTurbines(turb_coords, (turb_diam/2), 1)
hold off