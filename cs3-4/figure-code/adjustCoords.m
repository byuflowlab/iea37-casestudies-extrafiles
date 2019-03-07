% Nick Baker
% IEA37 cs3&4
% Adjust the original coordinates (Eastings and Northings) to zero out
clear, close all;

addpath('../support-files/')
row = 0;
col = 1; % Start on the second column (it can't read letters)

% Read .csv files
%[File3a, File3b, File4a, File4b, File4c] = readBorsseleBoundaries(row,col);
f1 = csvread('borsele-boundary-iiia(orig).csv',row,col);
% Adjust
%[File3a, File3b, File4a, File4b, File4c] = offsetBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c);
f1 = makeOffset(f1);
% Plot farm boundaries to check
%plotBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c);
% Given five boundary regions, plots them on the same graph
    hold on
    plotClosedBoundary(f1)
    % Zero out and scale our upper limits
    ylim([0 2.25e4])
    xlim([0 1.9e4])
    % Make sure it's proportional
    axis square
    hold off
% write to .csv files
%writeBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c);

%--- Changes the input coordinates by our Easting and Northing offset ---%
function [coordsOut] = makeOffset(coordsIn)
    coordOffset = [484178.50, 5716513.50];  % from pts 4 and 46

    coordsOut(:,1) = coordsIn(:,1) - coordOffset(1);    % Change x
    coordsOut(:,2) = coordsIn(:,2) - coordOffset(2);    % Change y
end

%--- Appends the first points in the list for the closed plot ---%
function [] = plotClosedBoundary(coords)
    plot([coords(:,1);coords(1,1)], [coords(:,2);coords(1,2)])
end

function [f1,f2,f3,f4,f5] = readBorsseleBoundaries(row,col)
    f1 = csvread('borsele-boundary-iiia(orig).csv',row,col);
    f2 = csvread('borsele-boundary-iiib(orig).csv',row,col);
    f3 = csvread('borsele-boundary-iva(orig).csv',row,col);
    f4 = csvread('borsele-boundary-ivb(orig).csv',row,col);
    f5 = csvread('borsele-boundary-ivc(orig).csv',row,col);
end

function [f1, f2, f3, f4, f5] = offsetBorsseleBoundaries(f1,f2,f3,f4,f5)
    f1 = makeOffset(f1);
    f2 = makeOffset(f2);
    f3 = makeOffset(f3);
    f4 = makeOffset(f4);
    f5 = makeOffset(f5);
end

function [] = plotBorsseleBoundaries(b1,b2,b3,b4,b5)
    % Given five boundary regions, plots them on the same graph
    hold on
    plotClosedBoundary(b1)
    plotClosedBoundary(b2)
    plotClosedBoundary(b3)
    plotClosedBoundary(b4)
    plotClosedBoundary(b5)
    % Zero out and scale our upper limits
    ylim([0 2.25e4])
    xlim([0 1.9e4])
    % Make sure it's proportional
    axis square
    hold off
end

function writeBorsseleBoundaries(File3a,File3b,File4a,File4b,File4c)
    dlmwrite('../support-files/borsele-boundary-iiia(adj).csv',File3a,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-iiib(adj).csv',File3b,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-iva(adj).csv',File4a,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-ivb(adj).csv',File4b,'precision','%.3f');
    dlmwrite('../support-files/borsele-boundary-ivc(adj).csv',File4c,'precision','%.3f');
end