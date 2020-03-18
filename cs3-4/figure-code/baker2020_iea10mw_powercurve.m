% Nicholas F. Baker
% IEA Task 37
% cs3&4
% NREL 10MW Power Curve figure
% 31 May 2018
% Updated 15 Jun (Adjusted for NREL 3.35 MW onshore reference turbine)
% Adjusted 27 Jan 2020 for IEA 10MW offshore reference turbine
clear, close all;

PlotColor = [0,0.4470,0.7410]; % Blue
%PlotColor = [0.6350, 0.0780, 0.1840]; % Red
%PlotColor = [0,0.4470,0.7410]; % Yellow

%{
% NREL 5MW offshore reference turbine data
nMaxWind = 15;
nMaxPwr = 5;
Vci = 3;        % (m/s) Cut-In Wind Speed
Vrws = 11.4;    % (m/s) Rated Wind Speed
Vco = 25;       % (m/s) Cut-Out Wind Speed
PwrRtng = 5;    % MW.
%}

%{
% NREL 3.35MW onshore reference turbine data
nMaxWind = 15;  % For figure (r) lateral limit 
nMaxPwr = 3.35;
Vci = 4;        % (m/s) Cut-In Wind Speed
Vrws = 9.8;     % (m/s) Rated Wind Speed
Vco = 25;       % (m/s) Cut-Out Wind Speed
PwrRtng = 3.35; % MW.
%}

% NREL 10MW onshore reference turbine data
nMaxWind = 30;  % For figure (r) lateral limit 
nMaxPwr = 10;
Vci = 4;        % (m/s) Cut-In Wind Speed
Vrws = 11;     % (m/s) Rated Wind Speed
Vco = 25;       % (m/s) Cut-Out Wind Speed
PwrRtng = 10; % MW.

xValues = linspace(0,nMaxWind, 500);
yValues = linspace(0,nMaxPwr*1.2, 500);

for i = 1:length(xValues)
    if (xValues(i) < Vci) % If we haven't hit the cut-in speed
        yValues(i) = 0;    % Give no power
    elseif (Vci <= xValues(i)) && (xValues(i) <= Vrws)
        yValues(i) = PwrRtng * ((xValues(i) - Vci) / (Vrws - Vci))^3;
    elseif (xValues(i) > Vco) % If we pass the cut-out speed
        yValues(i) = 0;    % Give no power
    else    % If we're faster than the full rated wind speed, use the full 
        yValues(i) = PwrRtng;
    end
end

%%%%%%%%% adjust below
% Plot the farm
figure('DefaultAxesFontSize',14);
p1 = plot(xValues,yValues, 'Color', PlotColor, 'LineWidth', 2); % Draw our curve
axis([0 nMaxWind 0 nMaxPwr*1.1]);                 % Make all plots the same size
%nNumTicks = fix((2*PlotDimen)/ (nMinSpaceDiam*TurbDiam))/2;
xticks(0:5:nMaxWind);
yticks(0:2:nMaxPwr);
ylabel('Power (MW)');
xlabel('Wind Speed (m/s)');
%title('Power Curve');
set(gca,'box','off')
grid off
%saveas(gcf,'iea37-10mw-pcurve.pdf')
%pbaspect([1 1 1])                                                  % Make the plots actually square