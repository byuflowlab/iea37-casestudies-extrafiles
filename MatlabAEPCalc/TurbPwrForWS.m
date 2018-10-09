function [turb_pwr] = TurbPwrForWS(wind_speed_eff, turb_ci, turb_co, rated_ws, rated_pwr)
    % Power curve of our used turbine, made a function to enable vectorization
    
    % If we're between the cut-in and rated wind speeds
    if ((turb_ci <= wind_speed_eff) && (wind_speed_eff < rated_ws))
        % Calculate the curve's power
        turb_pwr = rated_pwr * ((wind_speed_eff-turb_ci) / (rated_ws-turb_ci))^3;
    % If we're between the rated and cut-out wind speeds
    elseif ((rated_ws <= wind_speed_eff) && (wind_speed_eff < turb_co))
        % Produce the rated power
        turb_pwr = rated_pwr;
    end
end