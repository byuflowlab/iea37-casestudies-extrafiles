function [] = writeWindroseYAML(file_name, bin_degs, weib_vars, bin_freqs)
    %-- Function to write the turbine locations to a .yaml file in the iea37 format

    % Header info
    File.title = 'IEA Wind Task 37 Wind Plant Ontology version 0.1';
    File.description = 'Wind resource conditions for a IEA37 WFLO case studies 3 and 4.';
    File.definitions.wind_inflow.type = 'object';
    File.definitions.wind_inflow.description = 'inflow for current wind conditions';
    
    % Direction bins (in degrees)
    File.definitions.wind_inflow.properties.direction.id = 'wind_direction';
    File.definitions.wind_inflow.properties.direction.type = 'number';
    File.definitions.wind_inflow.properties.direction.description = 'The wind direction in degree, with North as the 0. 120 bins.';
    File.definitions.wind_inflow.properties.direction.units = 'deg';
    File.definitions.wind_inflow.properties.direction.bins = round(bin_degs,0);
    File.definitions.wind_inflow.properties.direction.minimum = 0.0;
    File.definitions.wind_inflow.properties.direction.maximum = 360.0;
    
    % Weibull parameters (Lambda, k) 
    File.definitions.wind_inflow.properties.speed.type = 'number';
    File.definitions.wind_inflow.properties.speed.desctiption = 'wind speeds for each bin, given as Weibull paramaters [Lambda, k]';
    File.definitions.wind_inflow.properties.speed.default = weib_vars;
    File.definitions.wind_inflow.properties.speed.units = 'm/s';

    % Turbulence intensity
    File.definitions.wind_inflow.properties.ti.type = 'object';
    File.definitions.wind_inflow.properties.ti.description = 'turbulence intensity';
    File.definitions.wind_inflow.properties.ti.default = 0.075;

    % Probability numbers
    File.definitions.wind_inflow.properties.probability.type = 'number';
    File.definitions.wind_inflow.properties.probability.description = 'Wind directional frequency distribution for 120 bins of wind rose';
    File.definitions.wind_inflow.properties.probability.default = bin_freqs;
    
    WriteYaml(file_name,File);        % Put it all in the file
end