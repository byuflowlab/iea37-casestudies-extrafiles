% Nicholas F. Baker
% IEA Task
% Optimization Only Case Study
% NREL 5MW Power Curve figure
% 31 May 2018
% Updated 15 Jun (Adjusted fr NREL 3.35 MW onshore reference turbine)
clear, close all;

%{
% NREL 5MW offshore reference turbine data
nMaxWind = 15;
nMaxPwr = 5;
Vci = 3;        % (m/s) Cut-In Wind Speed
Vrws = 11.4;    % (m/s) Rated Wind Speed
Vco = 25;       % (m/s) Cut-Out Wind Speed
PwrRtng = 5;    % MW.
%}

% NREL 3.35MW onshore reference turbine data
nMaxWind = 15;  % For figure (r) lateral limit 
nMaxPwr = 3.35;
Vci = 3;        % (m/s) Cut-In Wind Speed
Vrws = 9.8;     % (m/s) Rated Wind Speed
Vco = 25;       % (m/s) Cut-Out Wind Speed
PwrRtng = 3.35; % MW.
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

p1 = plot(xValues,yValues, 'Color', [0,0.4470,0.7410], 'LineWidth', 2); % Draw our curve

axis([0 nMaxWind 0 nMaxPwr*1.2]);                 % Make all plots the same size
%nNumTicks = fix((2*PlotDimen)/ (nMinSpaceDiam*TurbDiam))/2;
xticks(0:1:nMaxWind);
ylabel('Power (MW)');
xlabel('Wind Speed (m/s)');
%title('Power Curve');
grid on
%pbaspect([1 1 1])                                                  % Make the plots actually square