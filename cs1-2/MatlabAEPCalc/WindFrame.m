function [frame_coords] = WindFrame(turb_coords, wind_dir_deg)
    %-- Convert map coordinates to downwind/crosswind coordinates.

    % Convert from meteorological polar system (CW, 0 deg.=N)
    % to standard polar system (CCW, 0 deg.=W)
    % Shift so North comes "along" x-axis, from left to right.
    wind_dir_deg = 270 - wind_dir_deg;
    % Convert inflow wind direction from degrees to radians
    wind_dir_rad = deg2rad(wind_dir_deg);

    % Constants to use below
    cos_dir = cos(-wind_dir_rad);
    sin_dir = sin(-wind_dir_rad);
    % Convert to downwind(x) & crosswind(y) coordinates
    frame_coords.x = (turb_coords.x * cos_dir) - (turb_coords.y * sin_dir);
    frame_coords.y = (turb_coords.x * sin_dir) + (turb_coords.y * cos_dir);
end