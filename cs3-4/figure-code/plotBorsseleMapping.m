function [] = plotBorsseleMapping(fBndryPts, indexPts)
    % Given one boundary regions, map a grid
    NewPoints = zeros(length(indexPts),2);
    NewPoints(1,:) = fBndryPts(indexPts(1),:);
    NewPoints(2,:) = fBndryPts(indexPts(2),:);
    NewPoints(3,:) = fBndryPts(indexPts(3),:);
    NewPoints(4,:) = fBndryPts(indexPts(4),:);
    hold on
    scatter(NewPoints(:,1),NewPoints(:,2),'k', 'filled');
    hold off
end