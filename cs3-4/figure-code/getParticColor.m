% Nick Baker
% Broken into a function 02 Dec 18
% Givena a Participant's assigned number, pass back the color attributes.

function [PlotColor] = getParticColor(color_num)
    switch color_num
        case 1  % Red
            PlotColor = [0.6350, 0.0780, 0.1840];
        case 2  % Yellow
            PlotColor = [0.9290, 0.6940, 0.1250];
        case 3  % Purple
            PlotColor = [0.4940, 0.1840, 0.5560];
        case 4  % Green
            PlotColor = [0.4660, 0.6740, 0.1880];
        case 5  % Blue
            PlotColor = [0,0.4470,0.7410];
        case 6 % Orange
            PlotColor = [0.8500, 0.3250, 0.0980];
        case 7 % Puke yellow
            PlotColor = [0.75, 0.75, 0];
        case 8 % loyolagreen
            PlotColor = 1/255*[0,104,87];
        case 9 % Pink
            PlotColor = 1/255*[255, 204, 255];
        case 10 % loyolagray
            PlotColor = 1/255*[200,200,200];
        case 11 % Brown
            PlotColor = 1/255*[153, 102, 51];
    end
end