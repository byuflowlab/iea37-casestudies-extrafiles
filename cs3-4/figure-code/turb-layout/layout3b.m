function [TurbCoords] = layout3b(BndryCoords)
    % This is non-generic code that will specifically only work for the area 
    % IIIb layout
    nNumTurbines = 11;
    TurbCoords = zeros(nNumTurbines, 2); % x,y coords for all turbines
    
    %---- Map odd boundary to new boundary ----%
    %-- Find corner turbine (transcribe in tip) --%
    TurbCoords(1,:) = BndryCoords(8,:);     % First turbine goes in right tip
    TurbCoords(11,:) = BndryCoords(1,:);   % Corner

    %-- Print Debug --%
    %hold on
    %plotBorsseleBoundary(BndryCoords);
    %plotBorselleTurbines(TurbCoords, 99, 1);
    %hold off
    %-- End of Print Debug --%
   
    %-- Map along edges --%
    %- Right edge -%
    % Get x-dist, y-dist
    distRRowX = abs((BndryCoords(8,1) - BndryCoords(1,1))/ 5);
    distRRowY = abs((BndryCoords(8,2) - BndryCoords(1,2))/ 5);
    % Place R-side turbines in space
    rTurbIndx = [1,3,5,7,9,11];
    for i = 2:(length(rTurbIndx)-1)
        TurbCoords(rTurbIndx(i),1) = TurbCoords(rTurbIndx(i-1),1) - distRRowX;
        TurbCoords(rTurbIndx(i),2) = TurbCoords(rTurbIndx(i-1),2) - distRRowY;
    end

    %- Left edge -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(2,1) - BndryCoords(3,1))/ 6);
    distLRowY = abs((BndryCoords(2,2) - BndryCoords(3,2))/ 6);
    %- Place Pt2 -%
    TurbCoords(2,1) = BndryCoords(3,1) - (distLRowX/2);
    TurbCoords(2,2) = BndryCoords(3,2) - (distLRowY/2);
    %- Traverse and place -%
    leftPt = BndryCoords(2,:);
    rightPt = TurbCoords(2,:);
    distNeeded = abs(pdist([BndryCoords(2,:);BndryCoords(3,:)],'euclidean'))/5;
    % Loop through the leg, marking as we go
    lTurbIndx = [2,4,6,8,10];
    for i = 2:length(lTurbIndx)          % Lay each Turbine on the boundary
        TurbCoords(lTurbIndx(i),:) = findNewPtOnLine(rightPt,leftPt, distNeeded);   
        rightPt = TurbCoords(lTurbIndx(i),:);
    end
end