function [] = plotAEPs(nSubNums, fValues, nFarmSize)
    %nSubNums = flipud(nSubNums);
    %fValues = flipud(fValues);
    
    nNumSubs = length(nSubNums);
    nRanks = 1:nNumSubs;
    PlotColor = [];
    
    % get the plot colors
    for i = nRanks
        tempColor = getParticColor(nSubNums(i));
        PlotColor = [PlotColor; tempColor];
    end
    
    
    % Plot the farm
    figure(1);
    hold on
        p1 = scatter(nRanks,fValues,70, PlotColor, 'filled'); % A dumy thing for our legend
        b = num2str(nSubNums); c = cellstr(b);
        dx = 0.1;
        dy = 0.2; % displacement so the text does not overlay the data points
        text(nRanks+dx, fValues+dy, c, 'Fontsize', 13);
        ylabel('Percentage Increase from Example, %')
        xlabel('Submission Rank')
        axis([1 (nNumSubs-1) 0 18])
    hold off
end