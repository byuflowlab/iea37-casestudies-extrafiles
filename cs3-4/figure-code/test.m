% Nick Baker
% iea37 cs3&4
% test file to mess with code
clear, close all;
addpath('../startup-files/')
%addpath('turb-layout/')
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

% Read in YAML files
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex-opt4.yaml');
%[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex-opt3.yaml');
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_dir_freq, wind_speed, wind_speed_freq] = getWindroseYaml(fname_wr);
num_speed_bins = length(wind_speed);
[IIIa, IIIb, IVa, IVb, IVc] = getBorsBoundariesYaml('iea37-boundary-cs4.yaml');
%[IIIa] = getBorsBoundaryYaml('iea37-boundary-cs3.yaml');

%hold on
%plotBorsseleBoundaries(IIIa, IIIb, IVa, IVb, IVc)
%plotBorsseleBoundary(IIIa)
%plotBorselleTurbines(turb_coords, (turb_diam/2), 1)
%hold off

%Plot a Weibull slice
figure(1)
hold on
plot(wind_speed,wind_speed_freq(3,:))
scatter(wind_speed,wind_speed_freq(3,:),4,'o','MarkerFaceColor',[0.6350 0.0780 0.1840])
xlabel('Wind Speed (m/s)')
ylabel('Frequency')
ax = gca;
ax.FontSize = 12;
hold off

% Test of boundary files
% boundary = ReadYaml('iea37-boundary-cs4.yaml');    % Pull our .yaml file into a struct
% vert_list = cell2mat(boundary(1).boundaries(1).IIIa(:,:));
% loc_list = cell2mat(boundary(1).location(1).utm(:,:))

%plotWindRoseFreq(deg2rad(wind_dir), wind_dir_freq)
%[y] = calcWeibull(10, wind_speed(4,2), wind_speed(4,1))    % A test
% calc the AEP at the every . wind speed
%[AEP] = calcAEPcs3(turb_coords, wind_freq, wind_speed, num_speed_bins, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
