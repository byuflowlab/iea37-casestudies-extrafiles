function [pwrDir] = DirPower(turb_coords, wind_dir_deg, wind_speed, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
    %-- Return the power produced by all turbines, at each wind speed bin.
    frame_coords = WindFrame(turb_coords, wind_dir_deg);    % Shift coordinate frame of reference to downwind/crosswind
    loss = GaussianWake(frame_coords, turb_diam);           % Use the Simplified Bastankhah Gaussian wake model for wake deficits
    wind_speed_eff = wind_speed*(1-loss);                   % Effective windspeed is freestream multiplied by wake deficits

    % Check to see if turbine produces power for experienced wind speed
    turb_pwr = arrayfun(@(vector) TurbPwrForWs(vector, turb_ci, turb_co, rated_ws, rated_pwr), wind_speed_eff);
    % The power for this wind direction is the sum of that produced by each turbine
    pwrDir = sum(turb_pwr);
end