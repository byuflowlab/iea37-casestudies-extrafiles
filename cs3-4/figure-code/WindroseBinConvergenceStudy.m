% Nick Baker
% 8 Apr 19
% IEA37 windrose Convergence Study
clear, close all
addpath('../startup-files/')
addpath('../../../iea37-wflo-casestudies/cs3-4/')
addpath('../../cs1-2/MatlabAEPCalc')
addpath('turb-layout/')
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

% Filenames
TurbineLocFile = 'iea37-ex-opt4.yaml'; % Example turbine layout
%TurbineFile = '../startup-files/iea37-10mw.yaml';

% Get startup data
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML(TurbineLocFile);               % Get Turbine Location
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);    % Get Turbine data
[wind_dir, wind_freq, wind_speed] = getWindroseYaml(fname_wr);                      % Get original Wind Rose
wind_speed = [wind_speed(:,2), wind_speed(:,1)];                                    % Make it [k, Lambda]

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

% Variables
numOrigBins = length(wind_dir);     % number of original bins
littleMinBin = 11;
bigMinBin = 40;
littleBinStep = -1;
bigBinStep = -10;
maxWindSpeed = 25;                  % m/s
littleMaxBin = 40;
bigMaxBin = 120;
testLittleBins =  littleMaxBin:littleBinStep:littleMinBin;
testBigBins =  bigMaxBin:bigBinStep:bigMinBin;
numLittleSteps = length(testLittleBins);
numBigSteps = length(testBigBins);
littleAEP = zeros(numLittleSteps,1);
bigAEP = zeros(numBigSteps,1);
littlePercError = zeros(numLittleSteps,1);
bigPercError = zeros(numBigSteps,1);

%-- Do the benchmark of 360 degrees --%
Freqs = (oldBins(:,4)/100)';  % Make a decimal from a percentage
[~,g,~,tempFreqs] = extrapolateFrequencies(Freqs,60);
[f, newDirs,newWeibVars] = extrapolateWeibull([oldBins(:,3),oldBins(:,2)],60,maxWindSpeed);
newSpeeds = newWeibVars(:,2);
stdAEP = calcAEPcs3(turb_coords, tempFreqs, newWeibVars, 60, newDirs, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);

%-- Loop through major bins --%
% for i = 1:length(testBigBins)
%     % Get Wind Slices
%     [newDirs, newWeibVars] = getWeibSlices(f,wind_speed, testBigBins(i),maxWindSpeed);
% 	[newFreqs] = getFreqSlices(g, testBigBins(i));
%     newSpeeds = newWeibVars(:,2);
%     testBigBins(i)
% 	% Calculate AEP with wind rose
%     bigAEP(i) = calcAEPcs3(turb_coords, newFreqs, newWeibVars, 1, newDirs, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
%     bigPercError(i) = (bigAEP(i)/stdAEP)*100;
% end
% 
% %--  Loop through refined bins --%
% for i = 1:length(testLittleBins)
%     % Get Wind Slices
%     [newDirs, newWeibVars] = getWeibSlices(f,wind_speed, testLittleBins(i),maxWindSpeed);
% 	[newFreqs] = getFreqSlices(g, testLittleBins(i));
%     newSpeeds = newWeibVars(:,2);
%     testLittleBins(i)
% 	% Calculate AEP with wind rose
%     littleAEP(i) = calcAEPcs3(turb_coords, newFreqs,  newWeibVars, 1, newDirs, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
%     littlePercError(i) = (littleAEP(i)/stdAEP)*100;
% end

% Concatenate lists
% AEP = [stdAEP; bigAEP; littleAEP];
% PercError = [100; bigPercError; littlePercError];
% bins = [360,testBigBins,testLittleBins];
% PerDiff = abs(100-PercError);
% data = [bins', PerDiff];

% % Graph AEP vs bin number
% figure(1)
% hold on
% plot(bins,PerDiff)
% scatter(bins,PerDiff)
% hold off
% xlabel('Number of Bins')
% ylabel('Percent AEP from 360 deg')
% xlim([0 120])
% axis square
% grid on