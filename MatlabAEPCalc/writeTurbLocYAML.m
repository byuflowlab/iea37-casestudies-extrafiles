function [] = writeTurbLocYAML(file_name, turb_coords, fname_turb, fname_wr, binned_AEP)
    %-- Function to write the turbine locations to a .yaml file in the iea37 format

    % Header info
    File.input_format_version = 0;
    File.title = 'IEA Wind Task 37 Combined Case Study 9 Turbine Farm';
    File.description = 'Result file of Nick Baker using simplified Gaussian-Wake';

    % Turbine info
    File.definitions.wind_plant.type = 'object';
    File.definitions.wind_plant.description = 'specific plant design including turbine selection and placement';
    File.definitions.wind_plant.properties.layout.type = 'array';
    File.definitions.wind_plant.properties.layout.items{1} = [];
    File.definitions.wind_plant.properties.layout.items{1}.ref = '"#/definitions/position"';
    File.definitions.wind_plant.properties.layout.items{2} = fname_turb;

    % Turbine Locations
    File.definitions.position.type = 'array';
        num_coords = length(turb_coords);
        x_coords = turb_coords(1:(num_coords/2));
        y_coords = turb_coords((num_coords/2 +1):num_coords);
    File.definitions.position.items.xc = round(x_coords, 6);
    File.definitions.position.items.yc = round(y_coords, 6);
    File.definitions.position.additionalItems = false;
    File.definitions.position.description = 'an array of x-coordinates [x0, x1, ...] and y-coordinates [y0, y1, ...] of wind turbine positions in cartesian coordinates';
    File.definitions.position.units = 'm';

    % Wake model
    File.definitions.plant_energy.type = 'object';
    File.definitions.plant_energy.description = 'energy production from simplified Bastankhah Gaussian wake model';
    File.definitions.plant_energy.properties.wake_model_selection.type = 'algorithm';
    File.definitions.plant_energy.properties.wake_model_selection.description = 'wake model used to calculate AEP';
    File.definitions.plant_energy.properties.wake_model_selection.items{1}.ref = '"iea37-aepcalc.py"';

    % WindRose
    File.definitions.plant_energy.properties.wind_resource_selection.type = 'object';
    File.definitions.plant_energy.properties.wind_resource_selection.description = 'specific wind resource used to calculate AEP';
    File.definitions.plant_energy.properties.wind_resource_selection.properties.type = 'array';
    File.definitions.plant_energy.properties.wind_resource_selection.items{1}.ref = fname_wr;

    % AEP
    File.definitions.plant_energy.properties.annual_energy_production.type = 'number';
    File.definitions.plant_energy.properties.annual_energy_production.description = 'binned and total (default) annual energy production for a wind plant given a layout and binned wind rose';
    File.definitions.plant_energy.properties.annual_energy_production.binned = round(binned_AEP',6);
    File.definitions.plant_energy.properties.annual_energy_production.default = round(sum(binned_AEP),6);
    File.definitions.plant_energy.properties.annual_energy_production.units = 'MWh';

    WriteYaml(file_name,File);        % Put it all in the file
end