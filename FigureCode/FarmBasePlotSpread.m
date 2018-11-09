% Nicholas F. Baker
% IEA Task
% Optimization Only Case Study
% Wind Farm (Base Case) layout and depiction (plot)
% 29 May 2018
% Modified 08 June 2018 (Pushes outter ring to field boundary, evenly spaces inner rings
% Updated 15 Jun 2018 (Adjusted to use NREL 3.35 MW onshore reference turbine)
clear all, close all;

%PlotColor = [0,0.4470,0.7410]; % Blue
PlotColor = [0.6350, 0.0780, 0.1840]; % Red
%PlotColor = [0.9290, 0.6940, 0.1250]; % Yellow
%PlotColor = [1,0,0];

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
FarmSize = 3; % 1 = 16 turbines, 2 = 36 turbines, 3 = 64 turbines.
bWriteGridsToFile = false; % true = write the turbine coordinates to a file, false = don't write.
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

if(nNumTurbField >= 1)                  % If there's more than one turbine in the field
    nNumTurbLeft = nNumTurbField - 1;   % Keep a decrementing total of number of turbines left in the field. Here we take out the center one.
    bKeepGoing = true;                  % Make the flag to know we still haven't mapped all our turbines.
else
    bKeepGoing = false;                 % Make the flag to say we've mapped all (1) of the turbines.
end
nRingCntr = zeros(1);                   % Counter for which ring we're on.

% Loops through each concentric ring, fitting as many turbines in the ring as will fit
while(bKeepGoing)
    nRingCntr = nRingCntr + 1;  % Start the next ring
    % Make ring and circumference
    fCircRadius(nRingCntr) = nRingCntr * fMinSpaceDist;   % Make the radius the minimum distance between hub centers, times whatever ring we're on.
    fCircCircum(nRingCntr) = 2 * pi * fCircRadius(nRingCntr);

    % Get arclength between min turbine distance.
    c = fMinSpaceDist;  % Chord length between neighboring rotors
    h(nRingCntr) = fCircRadius(nRingCntr) - sqrt( fCircRadius(nRingCntr)^2 - ((c^2)/4) ); % Perpendicular height from center of chord to ring.
    S(nRingCntr) = asin(c/ (h(nRingCntr) + (c^2 / (4*h(nRingCntr))) )) * (h(nRingCntr) + (c^2 / (4*h(nRingCntr))) ); % arc length between neighboring turbines

    % Get number of turbines that fit in that ring
    nNumTurb(nRingCntr) = fix(fCircCircum(nRingCntr)/S(nRingCntr)); % Finds number of arcs that will fit on the circumference.
    nNumTurbLeft = nNumTurbLeft - nNumTurb(nRingCntr);
    
    % If we see the next ring is sparse (<4 turbines), expand current ring to accomodate
    if (nNumTurbLeft < 4) && (nNumTurbLeft > 0)
        % Expand current ring's circumference to accomodate number we need
        fCircCircum(nRingCntr) = S(nRingCntr) * (nNumTurb(nRingCntr) + nNumTurbLeft);
        % Find the radius from our new circumference
        fCircRadius(nRingCntr) = fCircCircum(nRingCntr)/ (2*pi);
        nNumTurbLeft = 0;   % Note that we've plotted all turbines now
        
        % Expand the inside rings as well
        fNewCircRadius = fCircRadius(nRingCntr)/nRingCntr; % To expand the inner rings as well
        for i = 1:(nRingCntr-1)
            % Adjust circumference of each inside ring
            fCircRadius(i) = i * fNewCircRadius;        % Make the radius the minimum distance between hub centers, times whatever ring we're on.
            fCircCircum(i) = 2 * pi * fCircRadius(i);
        end
    end

    % Increment to the next ring
    if(nNumTurbLeft <= 0)
       bKeepGoing = false; 
    end
end

%%{
%--- Beginning of modification to expand outisde ring to boundary, and fill all 
% Last ring of field boundary is field limit (FieldRad)
fCircRadius(nRingCntr) = fFieldRad;
fNewCircRadius = fCircRadius(nRingCntr)/nRingCntr; % To expand the inner rings as well
for i = 1:(nRingCntr-1)
    % Adjust circumference of each inside ring
    fCircRadius(i) = i * fNewCircRadius;        % Make the radius the minimum distance between hub centers, times whatever ring we're on.
    fCircCircum(i) = 2 * pi * fCircRadius(i);
end
%----- end of modification -----%
%}
% Take excess turbines out of last circle
nNumInnerTurbs = ones(1);   % Initialize to only one turbine in the center
for i = 1:(nRingCntr-1)     % Loop through all the rings besides the last one
   nNumInnerTurbs = nNumInnerTurbs + nNumTurb(i);
end
nNumTurb(nRingCntr) = nNumTurbField - nNumInnerTurbs;

% Assign grid coords for each turbine, first at center
fGridCoords = zeros(nNumTurbField,2); % (x-coord,y-coord). Array of grid coordinates in the field.
nCurrTurb = ones(1);
fTheta = zeros(nRingCntr,1);      % Central angle for each ring

for i = 1:nRingCntr             % For each ring of our circle
    nCurrTurb = nCurrTurb + 1;  % Running counter for the current turbine we're on
    fGridCoords(nCurrTurb, 1) = fCircRadius(i); % Force the first turbine in the ring to be on x-axis, at radius
    fGridCoords(nCurrTurb, 2) = 0; % y-coord of first one (on axis)
    fTheta(i) = 360 / nNumTurb(i);
    nCurrAng = 0;               % Running total for our angle, start on the x-axis in Q1.
    
    for j = 2:nNumTurb(i)       % For each turbine on this ring (note fGridCoords[1,~] = (Radius,0)
        nCurrTurb = nCurrTurb + 1;  % Running counter for the current turbine we're on
        nCurrAng = nCurrAng + deg2rad(fTheta(i));
        fGridCoords(nCurrTurb, 1) = fCircRadius(i) * cos(nCurrAng); % x-coord
        fGridCoords(nCurrTurb, 2) = fCircRadius(i) * sin(nCurrAng); % y-coord
    end
end

%FieldRad = fCircRadius( length(fCircRadius) ) + TurbRad;
%%{
%------- Ploting the farm -----%
figure('DefaultAxesFontSize',14);
hold on
p1 = scatter(0,0, 1, PlotColor, 'filled'); % A dumy thing for our legend
p2 = viscircles([0,0], fFieldRad,'LineStyle','--', 'Color', 'k');   % Plot the field boundary
for i = 1:nNumTurbField     % Print the turbines as filled in 'circles'
    rectangle('Position',[fGridCoords(i, 1)-TurbRad, fGridCoords(i, 2)-TurbRad, TurbDiam, TurbDiam], 'FaceColor', PlotColor, 'EdgeColor', [0,0.4470,0.7410], 'Curvature', [1 1])
end
hold off
axis([-PlotDimen PlotDimen -PlotDimen PlotDimen]);                 % Make all plots the same size
nNumTicks = fix((2*PlotDimen)/ (nMinSpaceDiam*TurbDiam))/2;
xticks(linspace(-PlotDimen,PlotDimen, nNumTicks+1));
xticklabels({'-30','-20','-10','0','10','20','30'});
yticks(linspace(-PlotDimen,PlotDimen, nNumTicks+1));
yticklabels({'-30','-20','-10','0','10','20','30'});
ylabel('Rotor Diameters (D)');
pbaspect([1 1 1])                                                  % Make the plots actually square
grid on
legend([p1 p2], 'Turbine', 'Field boundary')
legend boxoff
%{
%  Plots for Jared with no grid, axis, numbers, or titles
figure(2)
hold on
rectangle('Position',[-TurbRad, -TurbRad, TurbDiam, TurbDiam], 'FaceColor', PlotColor, 'EdgeColor', PlotColor, 'Curvature', [1 1])
p1 = scatter(0,0, 1, PlotColor, 'filled'); % A dumy thing for our legend
p2 = viscircles([0,0], fFieldRad,'LineStyle','--', 'Color', 'k');   % Plot the field boundary
for i = 1:nNumTurbField     % Print the turbines as filled in 'circles'
    rectangle('Position',[fGridCoords(i, 1)-TurbRad, fGridCoords(i, 2)-TurbRad, TurbDiam, TurbDiam], 'FaceColor', PlotColor, 'EdgeColor', PlotColor, 'Curvature', [1 1])
end
hold off

axis([-PlotDimen PlotDimen -PlotDimen PlotDimen]);                 % Make all plots the same size
grid off
ax = gca;
ax.Visible = 'off';
% The axes is removed from PDF file, too.
print(gcf,'-dpdf','-r300');
pbaspect([1 1 1])                                                  % Make the plots actually square
%}
%--- Code to write coordinates to a .csv file
%mIndexList = (0:(nNumTurbField-1))';
%mCoordMat = [mIndexList, fGridCoords];
%{
if bWriteGridsToFile
    switch(FarmSize)
        case 0
            dlmwrite('BaseCoords9.csv', fGridCoords, 'delimiter', ',', 'precision', '%10.4f'); % Write the coordinates to a CSV file
        case 1
            dlmwrite('BaseCoords16.csv', fGridCoords, 'delimiter', ',', 'precision', '%10.4f'); % Write the coordinates to a CSV file
        case 2
            dlmwrite('BaseCoords36.csv', fGridCoords, 'delimiter', ',', 'precision', '%10.4f'); % Write the coordinates to a CSV file
        case 3
            dlmwrite('BaseCoords64.csv', fGridCoords, 'delimiter', ',', 'precision', '%10.4f'); % Write the coordinates to a CSV file
        otherwise
            error('Variable "FarmSize" not initilized properly');
    end
end
%}