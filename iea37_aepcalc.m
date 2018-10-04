% Nick Baker
% Created for IEA Task 37 Wind Farm Layout Optimization Case Study 1
% Create 17 Sept 18

clear all, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

% Get turbine data from .yaml
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex16.yaml');
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml('iea37-335mw.yaml');
[wind_dir, wind_freq, wind_speed] = getWindRoseYaml('iea37-windrose.yaml');

%--- Optimization portion (fmincon)


%--- Read in .yaml files ---%
function [turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(file_name)
    %-- Function to pull the relevant data for the given turbine .yaml file
    
    TurbineStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcuts for the .yaml paths
    defs = TurbineStruct(1).definitions(1);
    op_props = defs.operating_mode(1).properties(1);
    turb_props = defs.wind_turbine_lookup(1).properties(1);
    rotor_props = defs.rotor(1).properties(1);
    % Pull the data we need
    turb_ci = op_props.cut_in_wind_speed(1).default(1);
    turb_co = op_props.cut_out_wind_speed(1).default(1);
    rated_ws = op_props.rated_wind_speed(1).default(1);
    rated_pwr = turb_props.power(1).maximum(1);
    turb_diam = (rotor_props.radius(1).default(1)) * 2;
end

function [wind_dir, wind_freq, wind_speed] = getWindRoseYaml(file_name)
    %-- Function to pull the relevant data for the given wind rose .yaml file
    
    WindRoseStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcuts for the .yaml paths
    props = WindRoseStruct(1).definitions(1).wind_inflow(1).properties(1);
    % Pull the data we need
    wind_dir = props.direction(1).bins(1);
    wind_freq = props.probability(1).default(1);
    wind_speed = props.speed(1).default(1);
end

function [turb_coords, fname_turb, fname_wr] = getTurbLocYAML(file_name)
    %-- Function to pull the turbine locations from the example .yaml file
    
    fname_turb = 'no such file';
    fname_wr = 'no such file';
    
    TurbLocStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcut to reused .yaml path
    defs = TurbLocStruct(1).definitions(1);
    coords = defs.position(1).items(1);
    turb_coords.x = coords.xc(1);
    turb_coords.y = coords.yc(1);
   % Read Turbine locations into an array
   
   ref_list_turbs = defs.wind_plant(1).properties(1).layout(1);
   ref_list_wr = defs.plant_energy(1).properties(1).wind_resource_selection(1).properties(1);
   
   % Get size of lists
   num_turb_files = numel(ref_list_turbs.items);
   num_wr_files = numel(ref_list_wr.items);
   
   % Iterate until first reference not lead by a '#'
   for i = 1:num_turb_files
        ref = ref_list_turbs.items(i);  % Rip the reference list
        ref = ref{1}(1).x0x24ref;       % Take string out of cell and structs
        if ~contains(ref, '#')          % Search for first file w/o '#'
            fname_turb = ref;           % Note filename for other .yaml files
        end
   end
   for i = 1:num_wr_files
      ref = ref_list_turbs.items(i);    % Rip the reference list
        ref = ref{1}(1).x0x24ref;       % Take string out of cell and structs
        if ~contains(ref, '#')          % Search for first file w/o '#'
            fname_wr = ref;             % Note filename for other .yaml files
        end
   end
end

%--- DirPower
function [pwrDir] = DirPower(turb_coords, wind_dir_deg, wind_speed, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
%-- Return the power produced by each turbine.
    num_turb = len(turb_coords);
    
    frame_coords = WindFrame(turb_coords, wind_dir_deg);    % Shift coordinate frame of reference to downwind/crosswind
    loss = GaussianWake(frame_coords, turb_diam);           % Use the Simplified Bastankhah Gaussian wake model for wake deficits
    wind_speed_eff = wind_speed*(1-loss);                   % Effective windspeed is freestream multiplied by wake deficits
    turb_pwr = zeros(num_turb,1);                           % By default, the turbine's power output is zero

    % Check to see if turbine produces power for experienced wind speed
    for n = 1:num_turb
        % If we're between the cut-in and rated wind speeds
        if ((turb_ci <= wind_speed_eff(n)) && (wind_speed_eff(n) < rated_ws))
            % Calculate the curve's power
            turb_pwr(n) = rated_pwr * ((wind_speed_eff(n)-turb_ci) / (rated_ws-turb_ci))^3;
        % If we're between the rated and cut-out wind speeds
        elseif ((rated_ws <= wind_speed_eff(n)) && (wind_speed_eff(n) < turb_co))
            % Produce the rated power
            turb_pwr(n) = rated_pwr;
        end
    end
    % Sum the power from all turbines for this direction
    pwrDir = sum(turb_pwr);
 end
 
%--- GaussianWake
function [loss] = GaussianWake(frame_coords, turb_diam)
    %-- Return each turbine's total loss due to wake from upstream turbines
    % Equations and values explained in <iea37-wakemodel.pdf>
    num_turb = len(frame_coords);

    % Constant thrust coefficient
    CT = 4.0*1./3.*(1.0-1./3.);
    % Constant, relating to a turbulence intensity of 0.075
    k = 0.0324555;
    % Array holding the wake deficit seen at each turbine
    loss = np.zeros(num_turb);

    for i = 1:num_turb                   % Looking at each turb (Primary)
        loss_array = zeros(num_turb,1);  % Calculate the loss from all others
        for j = 1:num_turb               % Looking at all other turbs (Target)
            x = frame_coords.x(i) - frame_coords.x(j);   % Calculate the x-dist
            y = frame_coords.y(i) - frame_coords.y(j);   % And the y-offset
            if (x > 0)                   % If Primary is downwind of the Target
                sigma = k*x + turb_diam/sqrt(8);         % Calculate the wake loss
                % Simplified Bastankhah Gaussian wake model
                exponent = -0.5 * (y/sigma)^2;
                radical = 1 - CT/(8*sigma^2 / turb_diam^2);
                loss_array(j) = (1.-np.sqrt(radical)) * exp(exponent);
            end
            % Note that if the Target is upstream, loss is defaulted to zero
        % Total wake losses from all upstream turbs, using sqrt of sum of sqrs
        loss(i) = sqrt(sum(loss_array^2));
        end
    end
end

%--- WindFrame
function [frame_coords] = WindFrame(turb_coords, wind_dir_deg)
    %-- Convert map coordinates to downwind/crosswind coordinates.

    % Convert from meteorological polar system (CW, 0 deg.=N)
    % to standard polar system (CCW, 0 deg.=W)
    % Shift so North comes "along" x-axis, from left to right.
    wind_dir_deg = 270. - wind_dir_deg;
    % Convert inflow wind direction from degrees to radians
    wind_dir_rad = deg2rad(wind_dir_deg);

    % Constants to use below
    cos_dir = cos(-wind_dir_rad);
    sin_dir = sin(-wind_dir_rad);
    % Convert to downwind(x) & crosswind(y) coordinates
    frame_coords.x = (turb_coords.x * cos_dir) - (turb_coords.y * sin_dir);
    frame_coords.y = (turb_coords.x * sin_dir) + (turb_coords.y * cos_dir);
end

%--- Target AEP function
function [pwrDir] = iea37_aepcalc(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
    % Calculate the wind farm AEP
    num_bins = len(wind_freq);  % Number of bins used for our windrose

    % Power produced by the wind farm from each wind direction
    pwr_produced = np.zeros(num_bins);
    % For each wind bin
    for i = 1:num_bins
        % Find the farm's power for the current direction
        pwr_produced(i) = DirPower(turb_coords, wind_dir(i), wind_speed, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
    end
    %  Convert power to AEP
    hrs_per_year = 365*24;
    AEP = hrs_per_year * (wind_freq * pwr_produced);
    AEP = AEP / 1e6;  % Convert to MWh
end