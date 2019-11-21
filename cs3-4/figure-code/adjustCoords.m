% Nick Baker
% IEA37 cs3&4
% Adjust the original coordinates (Eastings and Northings) to zero out
clear, close all;

addpath('../support-files/')
addpath('turb-layout/')
row = 0;
col = 1;            % Start on the second column (it can't read letters)

%-- Read .csv files --%
[File3a, File3b, File4a, File4b, File4c] = readBorsseleBoundaries(row,col);
%f1 = csvread('borsele-boundary-iiia(orig).csv',row,col);
%-- Adjust so left and lowermost boundary values are 'zero' --%
[bp3a,bp3b,bp4a,bp4b,bp4c] = zeroBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c);
%f1 = makeZero(f1);

hold on
plotBorsseleBoundaries(bp3a,bp3b,bp4a,bp4b,bp4c);
hold off

%--- Changes the input coordinates by our Easting and Northing offset ---%
function [coordsOut] = makeZero(coordsIn)
    coordOffset = [484178.50, 5716513.50];  % from pts 4 and 46

    coordsOut(:,1) = coordsIn(:,1) - coordOffset(1);    % Change x
    coordsOut(:,2) = coordsIn(:,2) - coordOffset(2);    % Change y
end
function [f1,f2,f3,f4,f5] = readBorsseleBoundaries(row,col)
    f1 = csvread('borsele-boundary-iiia(orig).csv',row,col);
    f2 = csvread('borsele-boundary-iiib(orig).csv',row,col);
    f3 = csvread('borsele-boundary-iva(orig).csv',row,col);
    f4 = csvread('borsele-boundary-ivb(orig).csv',row,col);
    f5 = csvread('borsele-boundary-ivc(orig).csv',row,col);
end
function [f1, f2, f3, f4, f5] = zeroBorsseleBoundaries(f1,f2,f3,f4,f5)
    f1 = makeZero(f1);
    f2 = makeZero(f2);
    f3 = makeZero(f3);
    f4 = makeZero(f4);
    f5 = makeZero(f5);
end
function writeBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c)
    dlmwrite('../support-files/borsele-boundary-iiia.csv',File3a,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-iiib.csv',File3b,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-iva.csv',File4a,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-ivb.csv',File4b,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-ivc.csv',File4c,'precision','%.3f');
end
