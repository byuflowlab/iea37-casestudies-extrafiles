function [AEP] = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
    % Calculate the wind farm AEP

    % Power produced by the wind farm from each wind direction
    % For each wind bin
    pwr_produced = arrayfun(@(vector) DirPower(turb_coords, vector, wind_speed, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr), wind_dir);
    
    %  Convert power to AEP
    hrs_per_year = 365*24;
    AEP = hrs_per_year * (wind_freq .* pwr_produced);
    AEP = AEP / 1e6;  % Convert to MWh
end