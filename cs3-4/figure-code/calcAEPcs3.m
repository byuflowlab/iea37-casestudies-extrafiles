function [AEP] = calcAEPcs3(turb_coords, dir_freq, weibVars, num_speed_bins, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
    % Calculate the wind farm AEP
    
    % Power produced by the wind farm from each wind direction
    % For each wind bin times the probability
    numDirBins = length(weibVars);
    nMaxWindSpeed = 25;
    hrs_per_year = 365*24;
    pwr_produced_ws = zeros(num_speed_bins, numDirBins);
    pwr_produced_wswf = zeros(1,numDirBins);
    %speed_array = zeros(num_speed_bins, numDirBins);
    % Power(direction, windspeed)
    % f(i) * w(i)
    % Sum across both
    
    for i=1:numDirBins
        %speed_array = binWeibDist(weibVars(i,:), nMaxWindSpeed, num_speed_bins);
        %test_array(i,:) = speed_array(:,2)';
%        for j = 1:num_speed_bins
            % Sum power produced by all turbines for this wind direction and speed
            pwr_produced_ws(i,:) = arrayfun(@(vector) DirPower(turb_coords, vector, weibVars(i,1), turb_diam, turb_ci, turb_co, rated_ws, rated_pwr), wind_dir);
            % Multiply by probability that this wind speed happens
%            pwr_produced_ws(j,:) = pwr_produced_ws(j,:) * speed_array(j,2);
%        end
%        sum(pwr_produced_ws)
%        pwr_produced_wswf(i,:) = sum(pwr_produced_ws);
%        pwr_produced_wswf(i,:) = pwr_produced_ws .* dir_freq;        % Multiply by frequency this direction happens
    end
    
    %  Convert power to AEP
    AEP = hrs_per_year .* (pwr_produced_ws .* dir_freq);
    AEP = sum(AEP) / 1e6;  % Convert to MWh
end