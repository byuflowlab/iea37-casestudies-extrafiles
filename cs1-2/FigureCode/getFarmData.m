% Nick Baker
% 02 Dec 18
% For a given participant and farm size, returns the relevant wind farm data

function [turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, fig_name, farm_rad, plot_dimen] = getFarmData(nParNum, nFarmSize)
    switch(nFarmSize)
        case 0
            farm_rad = 900;
            plot_dimen = 1200;
            switch(nParNum)
                case 0  % Example layout
                    fname_turb_loc = strcat('iea37-ex9.yaml');
                    fig_name = strcat('iea37-ex9.pdf');
                otherwise % All other participants
                    fname_turb_loc = strcat('iea37-par', num2str(nParNum), '-opt9.yaml');
                    fig_name = strcat('iea37-opt9-par', num2str(nParNum), '.pdf');
            end
        case 1
            farm_rad = 1300;
            plot_dimen = 1600;
            switch(nParNum)
                case 0  % Example layout
                    fname_turb_loc = strcat('iea37-ex16.yaml');
                    fig_name = strcat('iea37-ex16.pdf');
                otherwise % All other participants
                    fname_turb_loc = strcat('iea37-par', num2str(nParNum), '-opt16.yaml');
                    fig_name = strcat('iea37-opt16-par', num2str(nParNum), '.pdf');
            end
        case 2
            farm_rad = 2000;
            plot_dimen = 2500;
            switch(nParNum)
                case 0  % Example layout
                    fname_turb_loc = strcat('iea37-ex36.yaml');
                    fig_name = strcat('iea37-ex36.pdf');
                otherwise % All other participants
                    fname_turb_loc = strcat('iea37-par', num2str(nParNum), '-opt36.yaml');
                    fig_name = strcat('iea37-opt36-par', num2str(nParNum), '.pdf');
            end
        case 3
            farm_rad = 3000;
            plot_dimen = 3500;
            switch(nParNum)
                case 0  % Example layout
                    fname_turb_loc = strcat('iea37-ex64.yaml');
                    fig_name = strcat('iea37-ex64.pdf');
                otherwise % All other participants
                    fname_turb_loc = strcat('iea37-par', num2str(nParNum), '-opt64.yaml');
                    fig_name = strcat('iea37-opt64-par', num2str(nParNum), '.pdf');
            end
        otherwise
            error('Variable "FarmSize" not initilized properly');
    end

    % Get turbine data from .yaml
    [turb_coords, fname_turb, fname_wr] = getTurbLocYAML(fname_turb_loc);
    [turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
    [wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);
end