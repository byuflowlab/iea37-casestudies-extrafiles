function [newPt] = findNewPtOnLine(pt0,pt1,dist)
    % Given (x0,y0) and (x1,y1), finds (x2,y2) a distance dist from (x0,y0)
    % along line defined by (x0,y0) and (x1,y1).
    % https://math.stackexchange.com/questions/175896/finding-a-point-along-a-line-a-certain-distance-away-from-another-point

    distTot = abs(pdist([pt0;pt1],'euclidean'));
    t = dist/distTot;
    newPt(1) = ((1-t)*pt0(1)) + t*pt1(1);
    newPt(2) = ((1-t)*pt0(2)) + t*pt1(2);
end