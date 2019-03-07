function [] = plotSingleFarm(turb_coord, turb_diam, farm_rad, plot_dimen, color_num)
    % Plots a circular windfarm, given a farm radius, turbine diameter,
    % turbine coordinates (in .x and .y structured list) and a color

    turb_rad = turb_diam/2;
    num_turb = length(turb_coord.x);
    
    PlotColor = getParticColor(color_num);
    
    % Plot the farm
    figure(color_num);
    hold on
    %p1 = scatter(0,0, 1, PlotColor, 'filled'); % A dumy thing for our legend
    p2 = viscircles([0,0], farm_rad,'LineStyle',':', 'Color', 'k','LineWidth',.5);   % Plot the field boundary
    for i = 1:num_turb     % Print the turbines as filled in 'circles'
        rectangle('Position',[turb_coord.x(i)-turb_rad, turb_coord.y(i)-turb_rad, turb_diam, turb_diam], 'FaceColor', PlotColor, 'EdgeColor', PlotColor, 'Curvature', [1 1])
    end
    hold off
    axis([-plot_dimen plot_dimen -plot_dimen plot_dimen]);                 % Make all plots the same size
    %-- For grids and ticks: --%
    %{
        %nNumTicks = fix((2*plot_dimen)/ (nMinSpaceDiam*plot_dimen))/2;
        %xticks(linspace(-plot_dimen,plot_dimen, nNumTicks+1));
        %xticklabels({'-36','-24','-12','0','12','24','36'});
        %yticks(linspace(-plot_dimen,plot_dimen, nNumTicks+1));
        %yticklabels({'-36','-24','-12','0','12','24','36'});
        %ylabel('Rotor Diameters (D)');
        grid on
    %}
    %-- For no grids or ticks: --%
    %%{
        grid off;
        set(gca,'xtick',[])
        set(gca,'ytick',[])
        set(gca,'XColor','none')
        set(gca,'YColor','none')
        %set(gca,'Visible','off')
    %%}
    pbaspect([1 1 1])                                                  % Make the plots actually square
    %legend([p1 p2], 'Turbine', 'Field boundary')
end