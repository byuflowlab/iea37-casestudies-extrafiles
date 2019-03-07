% Nicholas F. Baker
% IEA Task
% Optimization Only Case Study
% Wind Farm Depiction (plot)
% 31 May 2018
clear, close all;

PlotColor = [0,0.4470,0.7410]; % Blue
%PlotColor = [0.6350, 0.0780, 0.1840]; % Red
%PlotColor = [0,0.4470,0.7410]; % Yellow

%{
% Turbine Data for NREL 5MW
TurbDiam = 126.4;           % Rotor Diameter, needed for constraints
TurbRad = TurbDiam/2;
%}
% Turbine Data for NREL 3.35MW
TurbDiam = 130.0;           % Rotor Diameter, needed for constraints
TurbRad = TurbDiam/2;
nMinSpaceDiam = 5;          % Minimum number of diameters in between turbines
fMinSpaceDist = (nMinSpaceDiam * TurbDiam);% + (2*TurbRad); % Minimum distance between rotors - the number of turbine diameters, + the radius of each turbine.
FarmSize = 0; % 1 = 16 turbines, 2 = 36 turbines, 3 = 64 turbines.
PlotDimen = 3900;

switch(FarmSize)
    case 0
        nNumTurbField = 9;
        fFieldRad = 900;
    case 1
        nNumTurbField = 16;                 % Overall number of turbines in the field.
        fFieldRad = 1300;
    case 2
        nNumTurbField = 36;                 % Overall number of turbines in the field.
        fFieldRad = 2000;
    case 3
        nNumTurbField = 64;                 % Overall number of turbines in the field.
        fFieldRad = 3000;
    otherwise
        error('Variable "FarmSize" not initilized properly');
end

%----- Plot the farm -----%
figure('DefaultAxesFontSize',14);

hold on
rectangle('Position',[-TurbRad, -TurbRad, TurbDiam, TurbDiam], 'FaceColor', PlotColor, 'EdgeColor', [0,0.4470,0.7410], 'Curvature', [1 1])
p1 = scatter(0,0, 1, PlotColor, 'filled'); % A dumy thing for our legend
p2 = viscircles([0,0], fFieldRad,'LineStyle','--', 'Color', 'k');   % Plot the field boundary

txtOrigin = '\leftarrow Reference turbine at (0,0)';
text(3*TurbRad,50,txtOrigin, 'FontSize',13)
hold off

axis([-PlotDimen PlotDimen -PlotDimen PlotDimen]);                 % Make all plots the same size
nNumTicks = fix((2*PlotDimen)/ (nMinSpaceDiam*TurbDiam))/2;
xticks(linspace(-PlotDimen,PlotDimen, nNumTicks+1));
xticklabels({'-36','-24','-12','0','12','24','36'});
yticks(linspace(-PlotDimen,PlotDimen, nNumTicks+1));
yticklabels({'-36','-24','-12','0','12','24','36'});
ylabel('NREL 5MW Rotor Diameters (D)');
grid on
pbaspect([1 1 1])                                                  % Make the plots actually square
legend([p1, p2],'NREL 5MW', 'Field boundary', 'Location', 'SE')
legend boxoff
%{
%  Plots for Jared with no grid, axis, numbers, or titles
figure(2)
hold on
rectangle('Position',[-TurbRad, -TurbRad, TurbDiam, TurbDiam], 'FaceColor', [0,0.4470,0.7410], 'EdgeColor', [0,0.4470,0.7410], 'Curvature', [1 1])
p1 = scatter(0,0, 1, [0,0.4470,0.7410], 'filled'); % A dumy thing for our legend
p2 = viscircles([0,0], fFieldRad,'LineStyle','--', 'Color', 'k');   % Plot the field boundary
hold off

axis([-PlotDimen PlotDimen -PlotDimen PlotDimen]);                 % Make all plots the same size
grid off
ax = gca;
ax.Visible = 'off';
pbaspect([1 1 1])                                                  % Make the plots actually square
%}