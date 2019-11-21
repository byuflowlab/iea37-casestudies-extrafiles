% Given coordiantes, calculates the interior area of a polygon
% Note: 1) shape must not be self intersecting.
%       2) coordinates must be clockwise or cc.
function [fArea] = calcPolygonArea(xCoords, yCoords)
    numCoords = length(xCoords);
    fSum = 0;
    
    % Sum from first to penultimate
    for i = 1:(numCoords-1)
        fFrstTrm = (xCoords(i) * yCoords(i+1));
        fScndTrm = (yCoords(i) * xCoords(i+1));
        fSum = fSum + (fFrstTrm - fScndTrm);
    end
    
    % Sum the last point (which wraps to the first)
    fFrstTrm = (xCoords(numCoords) * yCoords(1));
    fScndTrm = (yCoords(numCoords) * xCoords(1));
    fSum = fSum + (fFrstTrm - fScndTrm);
        
    % Div by 2 for answer
    fArea = abs(fSum / 2);
end