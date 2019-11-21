function [] = plotBorselleTurbines(TurbCoords, TurbRad, color)
    nNumTurbField = length(TurbCoords);
    TurbDiam = 2*TurbRad;
    for i = 1:nNumTurbField     % Print the turbines as filled in 'circles'
        rectangle('Position',[TurbCoords(i, 1)-TurbRad, TurbCoords(i, 2)-TurbRad, TurbDiam, TurbDiam], 'FaceColor', getParticColor(color), 'EdgeColor', getParticColor(color), 'Curvature', [1 1])
        %text(TurbCoords(i, 1),TurbCoords(i, 2),num2str(i),'FontSize',20,'Color','k');
    end
end