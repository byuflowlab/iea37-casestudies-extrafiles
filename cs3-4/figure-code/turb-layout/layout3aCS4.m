function [TurbCoords] = layout3aCS4(BndryCoords)
    % This is non-generic code that will specifically only work for the area 
    % IIIa layout
    nNumTurbines = 31;
    TurbCoords = zeros(nNumTurbines, 2); % x,y coords for all turbines
    
    %---- Map odd boundary to new boundary ----%
    %-- Find corner turbine (transcribe in tip) --%
    TurbCoords(1,:) = BndryCoords(1,:);     % First turbine goes in right tip
    TurbCoords(2,:) = BndryCoords(16,:);   % Corner
    TurbCoords(31,:) = BndryCoords(7,:);   % Corner
    TurbCoords(26,:) = BndryCoords(9,:);   % Corner
    TurbCoords(3,:) = BndryCoords(10,:);   % Corner

    %-- Print Debug --%
    %hold on
    %plotBorsseleBoundary(BndryCoords);
    %plotBorselleTurbines(TurbCoords, 99, 1);
    %hold off
    %-- End of Print Debug --%
    
    %-- Map along edges --%
    %- Right edge -%
%%%% USE IVA MAPPING!!!
    nNumLRows = 6;
    %- Left edge -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(9,1) - BndryCoords(10,1))/ (nNumLRows-1));
    distLRowY = abs((BndryCoords(9,2) - BndryCoords(10,2))/ (nNumLRows-1));
    % Place L-side turbines in space
    lTurbIndx = [3,7,11,16,21,26];
    for i = 2:(length(lTurbIndx)-1)
        TurbCoords(lTurbIndx(i),1) = TurbCoords(lTurbIndx(i-1),1) - distLRowX;
        TurbCoords(lTurbIndx(i),2) = TurbCoords(lTurbIndx(i-1),2) - distLRowY;
    end

    %- Right Edge -%
    distRight(1) = abs(pdist([BndryCoords(1,:);BndryCoords(2,:)],'euclidean'));
    distRight(2) = abs(pdist([BndryCoords(2,:);BndryCoords(3,:)],'euclidean'));
    distRight(3) = abs(pdist([BndryCoords(3,:);BndryCoords(4,:)],'euclidean'));
    distRight(4) = abs(pdist([BndryCoords(4,:);BndryCoords(5,:)],'euclidean'));
    distRight(5) = abs(pdist([BndryCoords(5,:);BndryCoords(6,:)],'euclidean'));
    distRight(6) = abs(pdist([BndryCoords(6,:);BndryCoords(7,:)],'euclidean'));
    distTotR = sum(distRight);
    rTurbIndx = [1,6,10,15,20,25,31];
    nNumRRows = length(rTurbIndx);
    distNeeded = distTotR / (nNumRRows-1); % Distance needed between each marker
    distRemain = distNeeded;                % Distance ramining till next marker
    legCntr = 1;                            % Start at first leg
    %-- Traverse and place --%
    leftPt = BndryCoords(1,:);
    rightPt = BndryCoords(2,:);
    distRemLeg = abs(pdist([leftPt;rightPt],'euclidean'));
    % Loop through all portions, marking as we go
    for ptCntr = 2:(nNumRRows-1)             % A turbine for each column
        while(distRemLeg < distRemain)          % If the remainder of our segment is smaller than we need
            distRemain = distRemain - distRemLeg;           % Take out that much distance
            legCntr = legCntr + 1;                          % Move to the next leg
            leftPt = BndryCoords(legCntr,:);                % Move our left point
            rightPt = BndryCoords(legCntr+1,:);
            distRemLeg = abs(pdist([leftPt;rightPt],'euclidean'));
        end
        TurbCoords(rTurbIndx(ptCntr),:) = findNewPtOnLine(leftPt,rightPt,distRemain);
        leftPt = TurbCoords(rTurbIndx(ptCntr),:);
        distRemLeg = distRemLeg - distRemain;       % Take out how much we used
        distRemain = distNeeded;                    % Reset howmuch we need
    end

    %- Bottom Edge -%
    distBottom(1) = abs(pdist([BndryCoords(7,:);BndryCoords(8,:)],'euclidean'));
    distBottom(2) = abs(pdist([BndryCoords(8,:);BndryCoords(9,:)],'euclidean'));
    distTotB = sum(distBottom);
    bTurbIndx = [31,30,29,28,27,26];
    nNumBRows = length(bTurbIndx);
    distNeeded = distTotB / (nNumBRows-1); % Distance needed between each marker
    distRemain = distNeeded;                % Distance ramining till next marker
    legCntr = 7;                            % Start at first leg
    %-- Traverse and place --%
    leftPt = BndryCoords(7,:);
    rightPt = BndryCoords(8,:);
    distRemLeg = abs(pdist([leftPt;rightPt],'euclidean'));
    % Loop through all portions, marking as we go
    for ptCntr = 2:(nNumBRows-1)             % A turbine for each column
        while(distRemLeg < distRemain)          % If the remainder of our segment is smaller than we need
            distRemain = distRemain - distRemLeg;           % Take out that much distance
            legCntr = legCntr + 1;                          % Move to the next leg
            leftPt = BndryCoords(legCntr,:);                % Move our left point
            rightPt = BndryCoords(legCntr+1,:);
            distRemLeg = abs(pdist([leftPt;rightPt],'euclidean'));
        end
        TurbCoords(bTurbIndx(ptCntr),:) = findNewPtOnLine(leftPt,rightPt,distRemain);
        leftPt = TurbCoords(bTurbIndx(ptCntr),:);
        distRemLeg = distRemLeg - distRemain;       % Take out how much we used
        distRemain = distNeeded;                    % Reset howmuch we need
    end

    %- Center turbines -%
    % Second Row %
    [TurbCoords] = fillRows([3,4,5,6], TurbCoords);
    % Third Row %
    [TurbCoords] = fillRows([7,8,9,10], TurbCoords);
    % Fourth Row %
    [TurbCoords] = fillRows([11,12,13,14,15], TurbCoords);
    % Fifth Row %
    [TurbCoords] = fillRows([16,17,18,19,20], TurbCoords);
    % Sixth Row %
    [TurbCoords] = fillRows([21,22,23,24,25], TurbCoords);
end
function [TurbCoords] = fillRows(indxList, TurbCoords)
    numTurbs = length(indxList);
    for i = 2:(numTurbs-1)
        % X-coord
        xFactor = (i-1)/(numTurbs-1);
        TurbCoords(indxList(i),1) = abs((TurbCoords(indxList(1),1) - TurbCoords(indxList(numTurbs),1)))*(xFactor) + TurbCoords(indxList(1),1);
        % Y-coord
        yFactor = 1-xFactor;
        TurbCoords(indxList(i),2) = abs((TurbCoords(indxList(1),2) - TurbCoords(indxList(numTurbs),2)))*(yFactor) + TurbCoords(indxList(numTurbs),2);
    end
end