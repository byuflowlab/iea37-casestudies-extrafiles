function [] = plotBorsseleSlices(fBndryPts)
    % Given one boundary regions, mark vertical slices at inflection points
    hold on
    % Delineate points
    for i=1:length(fBndryPts)
        line([fBndryPts(i,1) fBndryPts(i,1)], [fBndryPts(i,2) 0]);
    end
    % Draw vertical lines to zero (not above)
    
    hold off
end