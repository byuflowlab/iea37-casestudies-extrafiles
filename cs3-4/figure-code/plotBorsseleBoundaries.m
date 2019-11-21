function [] = plotBorsseleBoundaries(b1,b2,b3,b4,b5)
    % Given five boundary regions, plots them on the same graph
    %hold on
    plotClosedBoundary(b1,5,2, '-');
    plotClosedBoundary(b1,2,2, '--');
    plotClosedBoundary(b2,5,2, '-');
    plotClosedBoundary(b3,5,2, '-');
    plotClosedBoundary(b4,5,2, '-');
    plotClosedBoundary(b5,5,2, '-');
    % Zero out and scale our upper limits
    ylim([0 1.21e4])
    xlim([0 1.10e4])
    % Make sure it's proportional
    yticks(linspace(0,1.21e4,5));
    xticks(linspace(0,1.10e4,5));
    gca.YRuler.Exponent = 0;
    axis square
    axis off
    %hold off
end