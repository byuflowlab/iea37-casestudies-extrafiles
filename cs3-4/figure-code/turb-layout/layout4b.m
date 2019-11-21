% Assigns turbines to a generic layout igven a boundary

function [TurbCoords] = layout4b(BndryCoords)
    % This is non-generic code that will specifically only work for the area 
    % IVb layout
    nNumTurbine = 14;
    TurbCoords = zeros(nNumTurbine, 2); % x,y coords for all turbines
    % Make grid layout in triangle
    nNumRows = 6;
    
    %---- Map triangle boundary to new boundary ----%
    %-- Find top turbine (transcribe in tip) --%
    TurbCoords(1,:) = BndryCoords(1,:);     % First turbine goes in tip
    TurbCoords(14,:) = BndryCoords(2,:);   % Corner
    TurbCoords(11,:) = BndryCoords(3,:);   % Corner
   
    %-- Map along edges --%
    %- Right edge -%
    % Get x-dist, y-dist
    distRRowX = abs((BndryCoords(2,1) - BndryCoords(1,1))/ (nNumRows-1));
    distRRowY = abs((BndryCoords(2,2) - BndryCoords(1,2))/ (nNumRows-1));
    % Place R-side turbines in space
    rTurbIndx = [1,2,4,7,10,14];
    for i = 2:(nNumRows-1)
        TurbCoords(rTurbIndx(i),1) = TurbCoords(rTurbIndx(i-1),1) - distRRowX;
        TurbCoords(rTurbIndx(i),2) = TurbCoords(rTurbIndx(i-1),2) - distRRowY;
    end
    % Log weird 2nd turbine placement
    ScndR = TurbCoords(2,:);

    %-- Print Debug --%
    %hold on
    %plotBorsseleBoundary(BndryCoords);
    %plotBorselleTurbines(TurbCoords, 99, 1);
    %hold off
    %-- End of Print Debug --%

    %- Left edge -%
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(3,1) - BndryCoords(1,1))/ (nNumRows-1));
    distLRowY = abs((BndryCoords(3,2) - BndryCoords(1,2))/ (nNumRows-1));
    % Place R-side turbines in space
    lTurbIndx = [1,2,3,5,8,11];
    for i = 2:(nNumRows-1)
        TurbCoords(lTurbIndx(i),1) = TurbCoords(lTurbIndx(i-1),1) - distLRowX;
        TurbCoords(lTurbIndx(i),2) = TurbCoords(lTurbIndx(i-1),2) - distLRowY;
    end
    % Log weird 2nd turbine placement
    ScndL = TurbCoords(2,:);

    %- Bottom Row -%    
    % Get x-dist, y-dist
    distLRowX = abs((BndryCoords(3,1) - BndryCoords(2,1))/ 3);
    distLRowY = abs((BndryCoords(3,2) - BndryCoords(2,2))/ 3);
    % Place bottom turbines in space
    botTurbIndx = [11,12,13,14];
    for i = 2:(length(botTurbIndx)-1)
        TurbCoords(botTurbIndx(i),1) = TurbCoords(botTurbIndx(i-1),1) + distLRowX;
        TurbCoords(botTurbIndx(i),2) = TurbCoords(botTurbIndx(i-1),2) - distLRowY;
    end
    
    %- 2nd row to middle -%
    TurbCoords(2,1) = mean([ScndR(1);ScndL(1)]);
    TurbCoords(2,2) = mean([ScndR(2);ScndL(2)]);
    
    %- Move lower middle 2 -%
    TurbCoords(6,1) = mean([TurbCoords(5,1);TurbCoords(7,1)]);
    TurbCoords(6,2) = mean([TurbCoords(5,2);TurbCoords(7,2)]);
    TurbCoords(9,1) = mean([TurbCoords(8,1);TurbCoords(10,1)]);
    TurbCoords(9,2) = mean([TurbCoords(8,2);TurbCoords(10,2)]);
end