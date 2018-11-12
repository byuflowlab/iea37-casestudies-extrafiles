function [] = plotFarm(turb_coord, turb_diam, farm_rad, plot_dimen, color_num)
    % Plots a circular windfarm, given a farm radius, turbine diameter,
    % turbine coordinates (in .x and .y structured list) and a color

    turb_rad = turb_diam/2;
    num_turb = length(turb_coord.x);
    
    switch color_num
        case 1  % Red
            PlotColor = [0.6350, 0.0780, 0.1840];
        case 2  % Yellow
            PlotColor = [0.9290, 0.6940, 0.1250];
        case 3  % Purple
            PlotColor = [0.4940, 0.1840, 0.5560];
        case 4  % Green
            PlotColor = [0.4660, 0.6740, 0.1880];
        case 5  % Blue
            PlotColor = [0,0.4470,0.7410];
        case 6 % Orange
            PlotColor = [0.8500, 0.3250, 0.0980];
        case 7 % Puke yellow
            PlotColor = [0.75, 0.75, 0];
        case 8 % loyolagreen
            PlotColor = 1/255*[0,104,87];
        case 9 % Pink
            PlotColor = 1/255*[255, 204, 255];
        case 10 % loyolagray
            PlotColor = 1/255*[200,200,200];
        case 11 % Brown
            PlotColor = 1/255*[153, 102, 51];
    end
    
    % Plot the farm
    figure(1);
    hold on
    p1 = scatter(0,0, 1, PlotColor, 'filled'); % A dumy thing for our legend
    p2 = viscircles([0,0], farm_rad,'LineStyle','--', 'Color', 'k');   % Plot the field boundary
    for i = 1:num_turb     % Print the turbines as filled in 'circles'
        rectangle('Position',[turb_coord.x(i)-turb_rad, turb_coord.y(i)-turb_rad, turb_diam, turb_diam], 'FaceColor', PlotColor, 'EdgeColor', 'k', 'Curvature', [1 1])
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