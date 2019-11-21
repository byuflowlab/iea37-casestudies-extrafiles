% Nick Baker
% IEA37 cs4
% Places turbines in initial layout for cs4
clear, close all;

addpath('../support-files/')
addpath('turb-layout/')
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));
row = 0;
col = 0;            % Start on the second column (it can't read letters)
rotorRadius = 99;   % m, for 10 MW NREL turbin

%-- Read .csv files --%
[bp3a,bp3b,bp4a,bp4b,bp4c] = readBorsseleBoundaries(row,col);

%-- Shrink boundaires to accomodate rotor radii --%
[os3a,os3b,os4a,os4b,os4c] = offsetBorsseleBoundaries(bp3a,bp3b,bp4a,bp4b,bp4c,rotorRadius);

%-- Spread our Turbines --%
%[tc3a,tc3b,tc4a,tc4b,tc4c] = makeTurbineCoordsCS4(os3a,os3b,os4a,os4b,os4c);
%[tc3aCS3] = layout3aCS3(os3a);

%-- Plot farm boundaries to check --%
%hold on
%plotBorsseleBoundary(bp3a);
%plotBorselleTurbines(tc3aCS3, rotorRadius, 1);
%plotBorsseleBoundaries(bp3a,bp3b,bp4a,bp4b,bp4c);
%plotBorsseleBoundaries(os3a,os3b,os4a,os4b,os4c);
%plotBorselleTurbines(tc3a, rotorRadius, 1);
%plotBorselleTurbines(tc3b, rotorRadius, 1);
%plotBorselleTurbines(tc4a, rotorRadius, 1);
%plotBorselleTurbines(tc4b, rotorRadius, 1);
%plotBorselleTurbines(tc4c, rotorRadius, 1);
%hold off


%-- Present all turbines in .csv format --%
%allTurbLoc = round([tc3a;tc3b;tc4a;tc4b;tc4c],4);
%struct('position', allTurbLoc);
%WriteYaml('../support-files/cs4-baseline.yaml',allTurbLoc);
%WriteYaml('../support-files/cs3-baseline.yaml',round(tc3aCS3,4));
%WriteYaml('../support-files/cs3-outline.yaml',round(os3a,1));
WriteYaml('../support-files/cs3b-outline.yaml',round(os3b,1));
WriteYaml('../support-files/cs4a-outline.yaml',round(os4a,1));
WriteYaml('../support-files/cs4b-outline.yaml',round(os4b,1));
WriteYaml('../support-files/cs4c-outline.yaml',round(os4c,1));

function [f1,f2,f3,f4,f5] = readBorsseleBoundaries(row,col)
    f1 = csvread('../support-files/borsele-boundary-iiia.csv',row,col);
    f2 = csvread('../support-files/borsele-boundary-iiib.csv',row,col);
    f3 = csvread('../support-files/borsele-boundary-iva.csv',row,col);
    f4 = csvread('../support-files/borsele-boundary-ivb.csv',row,col);
    f5 = csvread('../support-files/borsele-boundary-ivc.csv',row,col);
end
function [tc3a,tc3b,tc4a,tc4b,tc4c] = makeTurbineCoordsCS4(os3a,os3b,os4a,os4b,os4c)
    [tc3a] = layout3aCS4(os3a);
    [tc3b] = layout3b(os3b);
    [tc4b] = layout4b(os4b);
    [tc4c] = layout4c(os4c);
    [tc4a] = layout4a(os4a);
end