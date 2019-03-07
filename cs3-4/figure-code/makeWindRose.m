% Nick Baker
% IEA37 cs3&4
% Wind Rose interpolation from discreet points
clear, close all

% Make matrix of points
Bins = zeros(13, 4);
%Bins(:,1) = linspace(0,360,13)'; % Makes 12 divisions of 0-360.
%Bins(end,:)=[];               % Delete last bin so it's 0-330  

% Hard code data from file
Bins(1,:) = [0,     8.65, 2.11, 5.1];
Bins(2,:) = [30.0,  8.86, 2.05, 4.3];
Bins(3,:) = [60.0,  8.15, 2.35, 4.3];
Bins(4,:) = [90.0,  9.98, 2.55, 6.6];
Bins(5,:) = [120.0, 11.35, 2.81, 8.9];
Bins(6,:) = [150.0, 10.96, 2.74, 6.5];
Bins(7,:) = [180.0, 11.28, 2.63, 8.7];
Bins(8,:) = [210.0, 11.50, 2.40, 11.5];
Bins(9,:) = [240.0, 11.08, 2.23, 12.0];
Bins(10,:) = [270.0, 10.94, 2.28, 11.1];
Bins(11,:) = [300.0, 11.27, 2.29, 11.4];
Bins(12,:) = [330.0, 10.55, 2.28, 9.6];
Bins(13,:) = [360.0,  8.65, 2.11, 5.1];      % Repeat of first point

% Convert to Polar (already in polar)
%[theta,rho] = cart2pol(Bins(:,1),Bins(:,2));
Bins(:,1) = deg2rad(Bins(:,1));

% Interpolate
pp = spline(Bins(:,1),Bins(:,3));           % Make the spline curve

% Plot
xPoints = linspace(0,2*pi)';                 % Do 100 points for a smooth curve
yPoints = ppval(pp,xPoints);                % Get specific Y-values
polarplot(xPoints,yPoints);                 % Plot it