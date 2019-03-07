function [] = plotOptHist(turb_coord, turb_diam, farm_rad, plot_dimen, color_num)
    % Plots a circular windfarm, given a farm radius, turbine diameter,
    % turbine coordinates (in .x and .y structured list) and a color

    turb_rad = turb_diam/2;
    num_turb = length(turb_coord.x);
    nMinSpaceDiam = 5;               % Minimum number of diameters in between turbines
    
    switch color_num
        case 0  % Blue
            PlotColor = [0,0.4470,0.7410];
        case 1  % Red
            PlotColor = [0.6350, 0.0780, 0.1840];
        case 2  % Yellow
            PlotColor = [0,0.4470,0.7410];
    end
    
    % Plot the farm
    figure(1);
    hold on
    p1 = scatter(0,0, 1, PlotColor, 'filled'); % A dumy thing for our legend
    p2 = viscircles([0,0], farm_rad,'LineStyle','--', 'Color', 'k');   % Plot the field boundary
    for i = 1:num_turb     % Print the turbines as filled in 'circles'
        rectangle('Position',[turb_coord.x(i)-turb_rad, turb_coord.y(i)-turb_rad, turb_diam, turb_diam], 'FaceColor', PlotColor, 'EdgeColor', [0,0.4470,0.7410], 'Curvature', [1 1])
    end
    hold off
    axis([-plot_dimen plot_dimen -plot_dimen plot_dimen]);                 % Make all plots the same size
    nNumTicks = fix((2*plot_dimen)/ (nMinSpaceDiam*plot_dimen))/2;
    xticks(linspace(-plot_dimen,plot_dimen, nNumTicks+1));
    xticklabels({'-36','-24','-12','0','12','24','36'});
    yticks(linspace(-plot_dimen,plot_dimen, nNumTicks+1));
    yticklabels({'-36','-24','-12','0','12','24','36'});
    ylabel('Rotor Diameters (D)');
    pbaspect([1 1 1])                                                  % Make the plots actually square
    grid on
    legend([p1 p2], 'Turbine', 'Field boundary')
end