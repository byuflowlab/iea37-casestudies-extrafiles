function [wind_dir, wind_dir_freq, wind_speed, wind_speed_freq] = getWindroseYaml(file_name)
    %-- Function to pull the relevant data for the given wind rose .yaml file (cs1&2)

    WindRoseStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcuts for the .yaml paths
    props = WindRoseStruct(1).definitions(1).wind_inflow(1).properties(1);
    % Pull the data we need
    wind_dir = cell2mat(props.direction(1).bins(:));
    wind_dir_freq = cell2mat(props.direction(1).frequency(:));
    wind_speed = cell2mat(props.speed(1).bins(:,:));
    wind_speed_freq = cell2mat(props.speed(1).frequency(:,:));
end