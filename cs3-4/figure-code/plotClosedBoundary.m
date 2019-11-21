%--- Appends the first points in the list for the closed plot ---%
function [] = plotClosedBoundary(coords, color, linewidth, style)
    plot([coords(:,1);coords(1,1)], [coords(:,2);coords(1,2)], 'Color', getParticColor(color), 'LineWidth', linewidth, 'LineStyle' ,style)
end