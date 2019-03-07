function [p1, parName] = plotAllFarm(turb_coord, turb_diam, farm_rad, plot_dimen, color_num, nNumPar)
    addpath(genpath('/Users/nbaker/Documents/MATLAB/subplot_tight/subplot_tight'));
    
    % Plots a subfigure circular windfarm, given a farm radius, turbine diameter,
    % turbine coordinates (in .x and .y structured list) and a color

    turb_rad = turb_diam/2;
    num_turb = length(turb_coord.x);
    
    PlotColor = getParticColor(color_num);
    
    % Plot the farm
    [m,n,p] = getSubfigPos(color_num, nNumPar);    % Get the subplot positioning
    subplot_tight(m,n,p);                                 % Plot this subplot
    hold on
    p1 = scatter(0,0, 1, PlotColor, 'filled', 'MarkerEdgeColor', 'k', 'Visible', 'off'); % A dumy thing for our legend
    p2 = viscircles([0,0], farm_rad,'LineStyle','--', 'Color', 'k','LineWidth',.2);   % Plot the field boundary
    for i = 1:num_turb     % Print the turbines as filled in 'circles'
        rectangle('Position',[turb_coord.x(i)-turb_rad, turb_coord.y(i)-turb_rad, turb_diam, turb_diam], 'FaceColor', PlotColor, 'EdgeColor', 'k', 'Curvature', [1 1])
    end
    hold off
    axis([-plot_dimen plot_dimen -plot_dimen plot_dimen]);                 % Make all plots the same size
    grid off;
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    
    pbaspect([1 1 1])                                                  % Make the plots actually square
    parName = strcat('sub', num2str(color_num), {' '});%'s', {' '} , num2str(num_turb),'-turbine layout');
    legend(p1, parName, 'location', 'southoutside')
    %p1.Visible = 'off';
    %legend('boxoff')           % Toggles the box around the legend
end