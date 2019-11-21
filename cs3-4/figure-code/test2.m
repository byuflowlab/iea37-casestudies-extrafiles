% Tests the turbine spacing for cs3
clear, close all;
addpath('../startup-files/')
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

fTurbRadius = 99;
shrink_factor = 1.75;
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex-opt3.yaml');
nNumTurbs = length(turb_coords);

bp3a = csvread('../support-files/borsele-boundary-iiia.csv',0,0);
bp3aTest = bp3a / shrink_factor;
turb_coords_test = turb_coords/shrink_factor;
plotBorsseleBoundary(bp3aTest);
plotBorselleTurbines(turb_coords_test, fTurbRadius, 1)


fArea(1) = calcPolygonArea(bp3aTest(:,1), bp3aTest(:,2));
fTurbArea = (pi*fTurbRadius^2)*nNumTurbs;

fDensity = (fTurbArea)/fArea;

distance = abs(pdist([turb_coords_test(22,:);turb_coords_test(17,:)],'euclidean'))/(fTurbRadius*2);