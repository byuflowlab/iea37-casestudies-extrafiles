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