function [wind_dir, wind_freq, wind_speed] = getWindRoseYaml(file_name)
    %-- Function to pull the relevant data for the given wind rose .yaml file
    
    WindRoseStruct = ReadYaml(file_name);    % Pull our .yaml file into a struct
    % Shortcuts for the .yaml paths
    props = WindRoseStruct(1).definitions(1).wind_inflow(1).properties(1);
    % Pull the data we need
    wind_dir = cell2mat(props.direction(1).bins(:));
    wind_freq = cell2mat(props.probability(1).default(:));
    wind_speed = props.speed(1).default(1);
end