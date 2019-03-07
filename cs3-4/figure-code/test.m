% Nick Baker
% iea37 cs3&4
% test file to mess with code
clear, close all;

printErikBounds();

function [] = printErikBounds()
    % Boundary points from Erik's .yaml files. Needed to visualize.
    xConst = [-0.316688,-0.965168,1.0,1.0,-0.316688];
    yConst = [1.0,-1.0,-0.144082,0.362028,1.0];
    xExc = [0.466179, 1.0, 0.466179];
    yExc = [1.0, -0.497298, 1.0];
    xBound = [0.521948, 0.676567, 0.200103, 0.429300, -0.870939, 0.521948];
    yBound = [-0.852978, 0.220158, 0.442276, 0.903162, 0.491392, -0.852978];

    hold on
    plot(xConst, yConst)
    plot(xExc, yExc)
    plot(xBound, yBound)
    axis square
    hold off
end