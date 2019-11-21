function [f,newDirsRad,newWeibVars] = extrapolateWeibull(oldWeibVars,numNewDirs,maxMagnitude)
    %--- Given a Weibull distribution of windspeeds around a windrose, splines and reslices new windspeed distributions.
    %-- <oldWeibVars> (:,1) k values (:,2) Lambda values or [k, Lambda] --%
    %-- <numNewDirs> number of new directions desired to slice --%
    
    %-- Setup --%
    %- Make times 3 for smooth splines on either side -%
    oldWeib3 = [oldWeibVars;oldWeibVars;oldWeibVars];
    %newWeibVars = zeros(numNewDirs,2);
    numOldDirs3 = length(oldWeibVars)*3;         % Multiply by three for our duplicates
    oldDirs3 = linspace((-2)*pi,4*pi-6*pi/numOldDirs3,numOldDirs3);
    newDirsRad = linspace(0,2*pi-(2*pi/numNewDirs),numNewDirs)';
    
    x = oldDirs3;
    y = linspace(0,maxMagnitude);                        % Y-axis is wind speeds (0 to ~25 m/s)
    z = zeros(length(x),length(y));                      % Z-axis is our Weibull distribution
    %figure(2)
    %hold on
    for i = 1:numOldDirs3
        z(i,:) = wblpdf(y,oldWeib3(i,1),oldWeib3(i,2));
    %    plot(y,z)
    end
    %hold off
    %-- Debug visualization --%
    %figure(1)
    %surf(x,y,z')
    %-- End debug visualization --%
    
    %- Spline for 3D surface -%
    [fpp]=csaps({x,y},z);          % Spline 3-D surface of all directions
    f = @(i,j) fnval(fpp,[i;j]);   % Turn the peicewise polynomial into a function
    
    [~,newWeibVars] = getWeibSlices(f,oldWeibVars, numNewDirs,maxMagnitude);
    %-- Debug visualization --%
    %xNew = linspace(0,2*pi-(2*pi/numNewDirs),numNewDirs)';
    %zNew = zeros(length(xNew),length(y));
    %for i = 1:numNewDirs
    %    zNew(i,:) = wblpdf(y,newWeibVars(i,1),newWeibVars(i,2));
    %end
    %figure(2)
    %surf(xNew,y,zNew')
end