% Nick Baker
% IEA37 cs3&4
% Shrink the boundary and turb locations by a factor so that the farm
% density is increased.
clear, close all;
addpath('../startup-files/')
addpath('../support-files/')
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

shrink_factor = 1.75;
% Read in Turbine Points
%[turb_coords3, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex-opt3-med.yaml');
[turb_coords4, fname_turb, ~] = getTurbLocYAML('iea37-ex-opt4-med.yaml');
[~, ~, ~, ~, turb_diam] = getTurbineYaml(fname_turb);
turb_rad = turb_diam/2;
% Shrink
%turb_coords3 = turb_coords3 / shrink_factor;
turb_coords4 = turb_coords4 / shrink_factor;
% Write Turb Points
%WriteYaml('opt3.yaml',round(turb_coords3,4));
%WriteYaml('opt4.yaml',round(turb_coords4,4));

% Read in boundary points
[File3a, File3b, File4a, File4b, File4c] = readBorsseleBoundaries(0,0);
[os3a,os3b,os4a,os4b,os4c] = offsetBorsseleBoundaries(File3a, File3b, File4a, File4b, File4c,99);
% Shrink
File3a = os3a/shrink_factor;
File3b = os3b/shrink_factor;
File4a = os4a/shrink_factor;
File4b = os4b/shrink_factor;
File4c = os4c/shrink_factor;
% Write Boundary Points
File.IIIa = round(File3a,1);
File.IIIb = round(File3b,1);
File.IVa = round(File4a,1);
File.IVb = round(File4b,1);
File.IVc = round(File4c,1);
%WriteYaml('small-boundary.yaml',File);

% Plot new farm with turbines
hold on
plotBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c);
%plotBorselleTurbines(turb_coords3, turb_rad, 1);
plotBorselleTurbines(turb_coords4, turb_rad, 1);
hold off


function [f1,f2,f3,f4,f5] = readBorsseleBoundaries(row,col)
    f1 = csvread('borsele-boundary-iiia.csv',row,col);
    f2 = csvread('borsele-boundary-iiib.csv',row,col);
    f3 = csvread('borsele-boundary-iva.csv',row,col);
    f4 = csvread('borsele-boundary-ivb.csv',row,col);
    f5 = csvread('borsele-boundary-ivc.csv',row,col);
end