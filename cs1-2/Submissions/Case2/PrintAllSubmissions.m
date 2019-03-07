% Nick Baker
% IEA37 Case Studies
% This code prints the participant proposed turbine locations in 2-d form.

clear, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/FigureCode'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/MatlabAEPCalc'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/Submissions/Case2'));
figuresdir = '/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/Figures/';

% Boolean values for functionality
bDispAEP = false;        % if True, shows calculations in command window
bWriteYaml = false;     % If True, writes calculated AEP data to a .yaml file
bDispFarm = true;      % If True, displays a fig of the wind farm
bSingleOrAll = false;    % If True, shows just one farm per fig. If false, displays all farms of a given size in a single fig with subfigures.
bSaveFig = true;       % If True, saves the WFL figure.
nNumPar = 5;           % Number of participants to graph
lgdColor = [];
lgdName = [];

%nFarmSize = 1;   % 0 = 9 turbines, 1 = 16 turbs, 2 = 36 turbs, 3 = 64 turbs
%nParmNum = 1;   % Participant number, 1-10.
for nFarmSize = 0:0          % Do all the farm sizes
    if (bDispFarm && ~(bSingleOrAll))    % If we are to display the farm, and it's got multiple figures
        figure('NumberTitle', 'off', 'Name', 'Case Study 2''s Layout Submissions');
        set(gcf,'Units','inches','Position',[0 0 8 11])
    end
    
    for nParNum = 1:10   % Do all the participants
        [turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, fig_name, farm_rad, plot_dimen] = getFarmData(nParNum, nFarmSize);
        nNumRtrs = length(turb_coords.x);   % Pulls the number of turbines by how many x-coordinates we have.

        % Get AEP data from turb locations
        if bDispAEP
            binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
            AEP = sum(binned_AEP)
        else
            binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
            AEP = sum(binned_AEP);
        end
        
        % Write the calculated data to a .yaml file
        if bWriteYaml       % If we need to write the farm AEP data to a file
            tc = [turb_coords.x,turb_coords.y];     % Strip the struct into a long array
            writeTurbLocYAML('TestAEP.yaml', tc, fname_turb, fname_wr, binned_AEP);
        end

        % Decide to display the actual farm or not.
        if bDispFarm
            if bSingleOrAll
                plotSingleFarm(turb_coords, turb_diam, farm_rad, plot_dimen, nParNum);
                % If we print a single farm, determine if we're saving the figure. 
                if bSaveFig
                    saveas(gcf,strcat(figuresdir, fig_name))            % Do it
                end
            else
                [tmpColor, tmpName] = plotAllFarm(turb_coords, turb_diam, farm_rad, plot_dimen, nParNum, nNumPar);
                % Add on the new stuff
                lgdColor = [lgdColor; tmpColor];
                lgdName = [lgdName; tmpName];
            end
        end
    end     % End of for-loop
    
    if (bDispFarm && ~(bSingleOrAll))    % If we are to display the farm, and it's got multiple figures
        plotSubfigLgd(lgdColor,lgdName,nNumPar);
        %sgtitle('Case Study 2''s Layout Submissions', 'FontWeight', 'bold')
        
        if bSaveFig                         % If we are to save the fig with subfigures
            fig_name = strcat('iea37-opt', num2str(nNumRtrs) ,'-', num2str(nNumPar) , 'plots.pdf');
            saveas(gcf,strcat(figuresdir, fig_name))            % Do it
        end
    end
end