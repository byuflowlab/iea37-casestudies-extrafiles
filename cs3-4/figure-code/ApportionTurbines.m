% Nick Baker
% IEA37 cs3&4
% Code to apportion turbines in each area based on their % of total area
% 18 Mar 2019
clear, close all;

addpath('../support-files/')
row = 0;
col = 0;
numBounds = 5;                  % Number of boundaries we're mapping
fArea = zeros(numBounds,1);     % Holds the actual area for each boundary
pArea = zeros(numBounds,1);     % Holds the percentage area for each boundary
numTurbs = zeros(numBounds,1);  % Holds the number of turbines in each boundary
switch(numBounds)
    case 5
        numTotTurbs = 81;
    case 1
        numTotTurbs = 25;
end

% Read .csv files
[File3a, File3b, File4a, File4b, File4c] = readBorsseleBoundaries(row,col);
%plotBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c);

% Find Area for each section
fArea(1) = calcPolygonArea(File3a(:,1), File3a(:,2));
fArea(2) = calcPolygonArea(File3b(:,1), File3b(:,2));
fArea(3) = calcPolygonArea(File4a(:,1), File4a(:,2));
fArea(4) = calcPolygonArea(File4b(:,1), File4b(:,2));
fArea(5) = calcPolygonArea(File4c(:,1), File4c(:,2));

% Find percent area for each section & correspond number of turbs per area
totArea = sum(fArea);   % Sum the areas
for i = 1:numBounds
    pArea(i) = fArea(i)/totArea;            % Find percentage area
    numTurbs(i) = numTotTurbs * pArea(i);   % Find percentage of turbines
end

numTurbs = round(numTurbs,0);               % Make # of turbs whole #s
numTurbs(2) = numTurbs(2) - 1;              % manual adjust for sum to match
sum(numTurbs)

function [f1,f2,f3,f4,f5] = readBorsseleBoundaries(row,col)
    f1 = csvread('borsele-boundary-iiia.csv',row,col);
    f2 = csvread('borsele-boundary-iiib.csv',row,col);
    f3 = csvread('borsele-boundary-iva.csv',row,col);
    f4 = csvread('borsele-boundary-ivb.csv',row,col);
    f5 = csvread('borsele-boundary-ivc.csv',row,col);
end
