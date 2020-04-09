% Nick Baker
% IEA37 cs3&4
% Wind Rose interpolation from discreet points
% 18 Mar 19
clear, close all
addpath('../support-files/')
%addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));
addpath(genpath('C:/Users/Captain Baker/Documents/MATLAB/YAMLMatlab/'));

% Make matrix of points
nNumOrigBins = 12;
nNumFinalBins = 360;
maxWindSpeed = 25;
oldBins = zeros(nNumOrigBins, 4); 

% Hard code data from original file
% [dirDeg, Lambda, k, freq%]
oldBins(1,:) =  [0,     8.65,  2.11, 5.1];
oldBins(2,:) =  [30.0,  8.86,  2.05, 4.3];
oldBins(3,:) =  [60.0,  8.15,  2.35, 4.3];
oldBins(4,:) =  [90.0,  9.98,  2.55, 6.6];
oldBins(5,:) =  [120.0, 11.35, 2.81, 8.9];
oldBins(6,:) =  [150.0, 10.96, 2.74, 6.5];
oldBins(7,:) =  [180.0, 11.28, 2.63, 8.7];
oldBins(8,:) =  [210.0, 11.50, 2.40, 11.5];
oldBins(9,:) =  [240.0, 11.08, 2.23, 12.0];
oldBins(10,:) = [270.0, 10.94, 2.28, 11.1];
oldBins(11,:) = [300.0, 11.27, 2.29, 11.4];
oldBins(12,:) = [330.0, 10.55, 2.28, 9.6];

%-- Extrapolate our frequencies --%
Freqs = (oldBins(:,4)/100)';  % Make a decimal from a percentage
[~,~,newDirs,newFreqs] = extrapolateFrequencies(Freqs,nNumFinalBins); % Extrapolate our points.

%-- Extrapolate our Weibull distributions --%
[~,~, newWeibVars] = extrapolateWeibull([oldBins(:,2),oldBins(:,3)],nNumFinalBins,maxWindSpeed);
newDirs = round(rad2deg(newDirs),0);
newWeibVars = round(newWeibVars,2);
%csvwrite('../support-files/new-WiebVars.csv',newWeibVars);
newFreqs = round(newFreqs,4)';
%csvwrite('../support-files/new-freqs.csv',newFreqs);
sum(newFreqs)
%-- Plot and write the data --%
% To plot the Wind Frequency distribution
%plotWindRoseFreq(newFreqs,newDirs);
%writeWindroseYAML('../support-files/new-iea37-cs3-windrose-test.yaml', newDirs, newWeibVars, newFreqs);