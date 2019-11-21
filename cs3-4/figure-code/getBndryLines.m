function [BndryEqArray] = getBndryLines(BndryPts)
    % Given a set of closed boundary points, returns the m,b values for all
    % connecting lines in y=mx+b format
    nNumPts = length(BndryPts);
    BndryEqArray = zeros(nNumPts,2);
    
    % Iterate through points and make lines
    for i = 1:(nNumPts-1)
        BndryEqArray(i,:) = polyfit([BndryPts(i,1), BndryPts(i+1,1)], [BndryPts(i,2), BndryPts(i+1,2)], 1);
    end
    
    % Do last one (wrapped around to first)
    BndryEqArray(nNumPts,:) = polyfit([BndryPts(nNumPts,1), BndryPts(1,1)], [BndryPts(nNumPts,2), BndryPts(1,2)], 1);
end