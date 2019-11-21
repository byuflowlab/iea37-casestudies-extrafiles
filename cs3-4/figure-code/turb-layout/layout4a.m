function [TurbCoords] = layout4a(BndryCoords)
    % This is non-generic code that will specifically only work for the area 
    % IVa layout
    nNumTurbines = 16;
    TurbCoords = zeros(nNumTurbines, 2); % x,y coords for all turbines
    % Make grid layout in triangle
    nNumRows = 4;
    nNumCols = 4;
    
    %---- Map odd boundary to new boundary ----%
    %-- Find top turbine (transcribe in tip) --%
    TurbCoords(1,:) = BndryCoords(1,:);     % First turbine goes in tip
    TurbCoords(4,:) = BndryCoords(4,:);   % Corner
    TurbCoords(16,:) = BndryCoords(5,:);   % Corner
    TurbCoords(13,:) = BndryCoords(6,:);   % Corner

    %-- Print Debug --%
    %hold on
    %plotBorsseleBoundary(BndryCoords);
    %plotBorselleTurbines(TurbCoords, 99, 1);
    %hold off
    %-- End of Print Debug --%
   
    %-- Map along edges --%
    %- Right edge -%
    % Get x-dist, y-dist
    distRRowX = abs((BndryCoords(5,1) - BndryCoords(4,1))/ (nNumRows-1));
    distRRowY = abs((BndryCoords(5,2) - BndryCoords(4,2))/ (nNumRows-1));
    % Place R-side turbines in space
    rTurbIndx = [4,8,12,16];
    for i = 2:(length(rTurbIndx)-1)
        TurbCoords(rTurbIndx(i),1) = TurbCoords(rTurbIndx(i-1),1) - distRRowX;
        TurbCoords(rTurbIndx(i),2) = TurbCoords(rTurbIndx(i-1),2) - distRRowY;
    end

    %- Left edge -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(6,1) - BndryCoords(1,1))/ (nNumRows-1));
    distLRowY = abs((BndryCoords(6,2) - BndryCoords(1,2))/ (nNumRows-1));
    % Place R-side turbines in space
    lTurbIndx = [1,5,9,13];
    for i = 2:(length(lTurbIndx)-1)
        TurbCoords(lTurbIndx(i),1) = TurbCoords(lTurbIndx(i-1),1) - distLRowX;
        TurbCoords(lTurbIndx(i),2) = TurbCoords(lTurbIndx(i-1),2) - distLRowY;
    end

    %- Bottom Row -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(6,1) - BndryCoords(5,1))/ (nNumCols-1));
    distLRowY = abs((BndryCoords(6,2) - BndryCoords(5,2))/ (nNumCols-1));
    % Place bottom turbines in space
    botTurbIndx = [13,14,15,16];
    for i = 2:(length(botTurbIndx)-1)
        TurbCoords(botTurbIndx(i),1) = TurbCoords(botTurbIndx(i-1),1) + distLRowX;
        TurbCoords(botTurbIndx(i),2) = TurbCoords(botTurbIndx(i-1),2) - distLRowY;
    end

    %- Top Edge -%
    distTop(1) = abs(pdist([BndryCoords(1,:);BndryCoords(2,:)],'euclidean'));
    distTop(2) = abs(pdist([BndryCoords(2,:);BndryCoords(3,:)],'euclidean'));
    distTop(3) = abs(pdist([BndryCoords(3,:);BndryCoords(4,:)],'euclidean'));
    distTopTot = sum(distTop);
    distNeeded = distTopTot / (nNumCols-1); % Distance needed between each marker
    distRemain = distNeeded;                % Distance ramining till next marker
    legCntr = 1;                            % Start at first leg
    %-- Traverse and place --%
    leftPt = BndryCoords(1,:);
    distRemLeg = abs(pdist([leftPt;BndryCoords(2,:)],'euclidean'));
    % Loop through all portions, marking as we go
    for ptCntr = 2:(nNumCols-1)             % A turbine for each column
        while(distRemLeg < distRemain)          % If the remainder of our segment is smaller than we need
            distRemain = distRemain - distRemLeg;           % Take out that much distance
            legCntr = legCntr + 1;                          % Move to the next leg
            leftPt = BndryCoords(legCntr,:);                % Move our left point
            distRemLeg = abs(pdist([leftPt;BndryCoords(legCntr+1,:)],'euclidean'));
        end
        TurbCoords(ptCntr,:) = findNewPtOnLine(leftPt,BndryCoords(legCntr+1,:),distRemain);
        distRemLeg = distRemLeg - distRemain;   % Take out how much we used
        distRemain = distNeeded;                % Reset howmuch we need
    end

    %-- Second Row --%
    TurbCoords(6,1) = abs((TurbCoords(2,1) - TurbCoords(14,1)))*(2/3) + TurbCoords(14,1);
    TurbCoords(6,2) = abs((TurbCoords(5,2) - TurbCoords(8,2)))*(2/3) + TurbCoords(8,2);
    TurbCoords(7,1) = abs((TurbCoords(3,1) - TurbCoords(15,1)))*(2/3) + TurbCoords(15,1);
    TurbCoords(7,2) = abs((TurbCoords(5,2) - TurbCoords(8,2)))*(1/3) + TurbCoords(8,2);

    %-- Third Row --%
    TurbCoords(10,1) = abs((TurbCoords(2,1) - TurbCoords(14,1)))*(1/3) + TurbCoords(14,1);
    TurbCoords(10,2) = abs((TurbCoords(5,2) - TurbCoords(8,2)))*(1/3) + TurbCoords(8,2);
    TurbCoords(11,1) = abs((TurbCoords(3,1) - TurbCoords(15,1)))*(1/3) + TurbCoords(15,1);
    TurbCoords(11,2) = abs((TurbCoords(9,2) - TurbCoords(12,2)))*(1/3) + TurbCoords(12,2);
end