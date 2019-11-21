% Nick Baker
% IEA37 Case Studies
% This code ranks AEP outputs and prints the percent difference from the example layout.

clear, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/FigureCode'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/MatlabAEPCalc'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/Submissions/Case2'));
figuresdir = '/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/Figures/';

% Boolean values for functionality
bDispAEP = false;           % if True, shows calculations in command window

% Constants for calculations
nNumPar = 12;               % Number of Particpanits we're comparing
AEP = zeros((nNumPar+1),1);% Add one for the example layout
nFarmSize = 3;              % 0 = 9 turbines, 1 = 16 turbs, 2 = 36 turbs, 3 = 64 turbs
%nParNum = 1;               % Participant number, 1-10.
for nParNum = 1:nNumPar     % Do all the participants
    % Get the turbine location and data for the participant's layout
    [turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, ~, ~, ~] = getFarmData(nParNum, nFarmSize);

    % Get AEP data from turb locations
    if bDispAEP
        binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
        AEP(nParNum) = sum(binned_AEP)
    else
        binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
        AEP(nParNum) = sum(binned_AEP);
    end
end     % End of for-loop

% Get example layout data, put it in the last spot in our array
[turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, fig_name, farm_rad, plot_dimen] = getFarmData(0, nFarmSize);
nNumRtrs = length(turb_coords.x);   % Pulls the number of turbines by how many x-coordinates we have.
binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
AEP(nParNum+1) = sum(binned_AEP);

% Sort highest to lowest
[sortedAEP,rankOrder] = sort(AEP,'descend');
%rankOrder = rankOrder';     % Order of participant AEP ranks

% Calculate and print differences in order
nExIndx = nNumPar+1;
fPerDiffEx = zeros(nNumPar,1);   % Array holding the percent difference from the example layout
for i = 1:nNumPar               % Do all the participants
    fDiff = sortedAEP(i) - AEP(nExIndx);          % Participant AEP from ex. AEP
    fPerDiffEx(i) = fDiff / AEP(nExIndx) * 100;   % Divide the difference by ex. AEP
end

%plotAEPs(rankOrder, sortedAEP)
fPerDiffEx = [fPerDiffEx;0]; % Add the example (0% diff)
plotAEPs(rankOrder, fPerDiffEx, nFarmSize)
%plotAEPs(1:nNumPar+1, fPerDiffEx, nFarmSize)

% Calculate percentage difference from highest, including example
nExIndx = nNumPar+1;
fPerDiffHigh = zeros(nNumPar,1);   % Array holding the percent difference from the example layout
for i = 1:(nNumPar+1)              % Do all the participants
    fDiff = sortedAEP(i) - sortedAEP(1);            % Participant AEP from highest AEP
    fPerDiffHigh(i) = fDiff / sortedAEP(1) * 100;   % Divide the difference by highest AEP
end