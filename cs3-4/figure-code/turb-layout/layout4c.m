function [TurbCoords] = layout4c(BndryCoords)
    % This is non-generic code that will specifically only work for the area 
    % IVc layout
    nNumTurbines = 9;
    TurbCoords = zeros(nNumTurbines, 2); % x,y coords for all turbines
    % Make grid layout in triangle
    nNumRows = 4;
    
    %---- Map triangle boundary to new boundary ----%
    %-- Find top turbine (transcribe in tip) --%
    TurbCoords(1,:) = BndryCoords(1,:);     % First turbine goes in tip
    TurbCoords(9,:) = BndryCoords(2,:);   % Corner
    TurbCoords(5,:) = BndryCoords(5,:);   % Corner
   
    %-- Map along edges --%
    %- Right edge -%
    % Get x-dist, y-dist
    distRRowX = abs((BndryCoords(2,1) - BndryCoords(1,1))/ 2);
    distRRowY = abs((BndryCoords(2,2) - BndryCoords(1,2))/ 2);
    % Place R-side turbines in space
    rTurbIndx = [1,6,9];
    for i = 2:(length(rTurbIndx)-1)
        TurbCoords(rTurbIndx(i),1) = TurbCoords(rTurbIndx(i-1),1) - distRRowX;
        TurbCoords(rTurbIndx(i),2) = TurbCoords(rTurbIndx(i-1),2) - distRRowY;
    end

    %- Left edge -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(5,1) - BndryCoords(1,1))/ 4);
    distLRowY = abs((BndryCoords(5,2) - BndryCoords(1,2))/ 4);
    % Place R-side turbines in space
    lTurbIndx = [1,2,3,4,5];
    for i = 2:(length(lTurbIndx)-1)
        TurbCoords(lTurbIndx(i),1) = TurbCoords(lTurbIndx(i-1),1) - distLRowX;
        TurbCoords(lTurbIndx(i),2) = TurbCoords(lTurbIndx(i-1),2) - distLRowY;
    end

    %- Bottom Row -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(5,1) - BndryCoords(2,1))/ 2);
    distLRowY = abs((BndryCoords(5,2) - BndryCoords(2,2))/ 2);
    % Place bottom turbines in space
    botTurbIndx = [5,8,9];
    for i = 2:(length(botTurbIndx)-1)
        TurbCoords(botTurbIndx(i),1) = TurbCoords(botTurbIndx(i-1),1) + distLRowX;
        TurbCoords(botTurbIndx(i),2) = TurbCoords(botTurbIndx(i-1),2) - distLRowY;
    end

    % One in middle
    TurbCoords(7,1) = abs(mean([TurbCoords(3,1);TurbCoords(9,1)]));
    TurbCoords(7,2) = abs(mean([TurbCoords(3,2);TurbCoords(9,2)]));

    %- Adjust bottom one -%
    % Get distance along bottom line
    distBottom(1) = abs(pdist([BndryCoords(2,:);BndryCoords(3,:)],'euclidean'));
    distBottom(2) = abs(pdist([BndryCoords(3,:);BndryCoords(4,:)],'euclidean'));
    distBottom(3) = abs(pdist([BndryCoords(4,:);BndryCoords(5,:)],'euclidean'));
    distBottomTot = sum(distBottom);
    distNeeded = distBottomTot / 2;    % divide by correct number
    distRemain = distNeeded;
    %-- Traverse and place --%
    % Subtract first portion
    distRemain = distRemain - distBottom(1);
    % Subtract second portion
    TurbCoords(8,:) = findNewPtOnLine(BndryCoords(3,:),BndryCoords(4,:),distRemain);
end