
clear, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));
addpath(genpath('/Users/nbaker/Documents/GitHub/iea37-casestudies-extrafiles/FigureCode'));

farmSize = 0;   % 0 = 9 turbines, 1 = 16 turbs, 2 = 36 turbs, 3 = 64 turbs
%participant_number = 9;     % Participant number, 1-10.
for participant_number = 5:5   % Do all the participants
    switch(farmSize)
        case 0
            fname_turb_loc = strcat('iea37-par', num2str(participant_number), '-opt9.yaml');
            fig_name = strcat('iea37-opt9-par', num2str(participant_number), '.pdf');
            farm_rad = 900;
            plot_dimen = 1200;
        case 1
            fname_turb_loc = strcat('iea37-par', num2str(participant_number), '-opt16.yaml');
            fig_name = strcat('iea37-opt16-par', num2str(participant_number), '.pdf');
            farm_rad = 1300;
            plot_dimen = 1600;
        case 2
            fname_turb_loc = strcat('iea37-par', num2str(participant_number), '-opt36.yaml');
            fig_name = strcat('iea37-opt36-par', num2str(participant_number), '.pdf');
            farm_rad = 2000;
            plot_dimen = 2500;
        case 3
            fname_turb_loc = strcat('iea37-par', num2str(participant_number), '-opt64.yaml');
            fig_name = strcat('iea37-opt64-par', num2str(participant_number), '.pdf');
            farm_rad = 3000;
            plot_dimen = 3500;
        otherwise
            error('Variable "FarmSize" not initilized properly');
    end

    % Get turbine data from .yaml
    [turb_coords, fname_turb, fname_wr] = getTurbLocYAML(fname_turb_loc);
    [turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
    [wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);
    nNumRtrs = length(turb_coords.x);   % Pulls the number of turbines by how many x-coordinates we have.

    % Get AEP data from turb locations
    binned_AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
    AEP = sum(binned_AEP)
    
    %tc = [turb_coords.x,turb_coords.y];     % Strip the struct into a long array
    %write_fname = strcat('iea37-par5-cc', num2str(participant_number), '.yaml');     % Save the cross comparison name
    %writeTurbLocYAML(write_fname, tc, fname_turb, fname_wr, binned_AEP);    % Write the computed AEP data

    %color_num = 2;  % 0 = blue, 1 = red, 2 = yellow, 3 = purple, 4 = green
    %clf     % Clear the figure for the next one
    %plotFarm(turb_coords, turb_diam, farm_rad, plot_dimen, participant_number)
    %saveas(gcf,fig_name)
end