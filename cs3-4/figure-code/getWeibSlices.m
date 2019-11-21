function [newDirsRad, newWeibVars] = getWeibSlices(f,oldWeibVars, numNewDirs,maxMagnitude)
    newWeibVars = zeros(numNewDirs,2);
    y = linspace(0,maxMagnitude);                        % Y-axis is wind speeds (0 to ~25 m/s)
    newDirsRad = linspace(0,2*pi-(2*pi/numNewDirs),numNewDirs)';
    
    %- Re-slice to new number -%
    xNew = linspace(0,2*pi-(2*pi/numNewDirs),numNewDirs)';
    zNew = zeros(length(xNew),length(y));
    for i = 1:numNewDirs
        for j = 1:length(y)
            zNew(i,j) = f(xNew(i),y(j));
        end
    end

    %- Extrapolate Weibull variables -%
    for i = 1:numNewDirs
        wblObj = fit(y',zNew(i,:)','(b/a)*(x/a)^(b-1)*exp(-(x/a)^b)', 'StartPoint', [oldWeibVars(1), oldWeibVars(2)]);
        newWeibVars(i,:) = [wblObj.a,wblObj.b];
    end
    %figure(1)
    %hold on
    %plot(y,zNew)
    %hold off
end