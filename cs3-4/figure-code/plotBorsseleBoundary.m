function [] = plotBorsseleBoundary(fBndryPts)
    % Given one boundary regions, plot it on the same scale as all 5
    hold on
    plotClosedBoundary(fBndryPts,5,1, '-');
    plotClosedBoundary(fBndryPts,2,2, '--');
    %labels = {'Point 1','Point 2','Point 3'};
    %plot(b1(:,1),b1(:,2),'o')
    %text(b1(:,1),b1(:,2),labels,'VerticalAlignment','bottom','HorizontalAlignment','right')
    % Zero out and scale our upper limits
    ylim([min(fBndryPts(:,2)) max(fBndryPts(:,2))])
    xlim([min(fBndryPts(:,1)) max(fBndryPts(:,1))])
    % Make sure it's proportional
    yticks(linspace(0,2.25e4,5));
    xticks(linspace(0,1.9e4,5));
    gca.YRuler.Exponent = 0;
    axis square
    axis off
    hold off
end