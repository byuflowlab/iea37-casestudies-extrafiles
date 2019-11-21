function [BndryPts1,BndryPts2,BndryPts3,BndryPts4,BndryPts5] = offsetBorsseleBoundaries(BndryPts1,BndryPts2,BndryPts3,BndryPts4,BndryPts5, fRotorRadius)
    %-- Shrinks boundaires to accomodate rotor radii --%
    nNumBounds = 5;  % How many boundaries we are working with
    nNumPoints = zeros(nNumBounds,1);
    nNumPoints(1) = length(BndryPts1);
    nNumPoints(2) = length(BndryPts2);
    nNumPoints(3) = length(BndryPts3);
    nNumPoints(4) = length(BndryPts4);
    nNumPoints(5) = length(BndryPts5);
    
    %-- Get the original boundary lines --%
    BndryEqArray1 = getBndryLines(BndryPts1);
    BndryEqArray2 = getBndryLines(BndryPts2);
    BndryEqArray3 = getBndryLines(BndryPts3);
    BndryEqArray4 = getBndryLines(BndryPts4);
    BndryEqArray5 = getBndryLines(BndryPts5);
    
    %-- Adjust towards the inside --%
    BndryPts1 = getInteriorLines(BndryEqArray1, BndryPts1, fRotorRadius);
    BndryPts2 = getInteriorLines(BndryEqArray2, BndryPts2, fRotorRadius);
    BndryPts3 = getInteriorLines(BndryEqArray3, BndryPts3, fRotorRadius);
    BndryPts4 = getInteriorLines(BndryEqArray4, BndryPts4, fRotorRadius);
    BndryPts5 = getInteriorLines(BndryEqArray5, BndryPts5, fRotorRadius);
end