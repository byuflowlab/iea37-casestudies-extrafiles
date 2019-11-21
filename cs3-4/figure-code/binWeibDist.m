% Nick Baker
% IEA37 cs3
% Written to check that the weibull curves sum correctly, even when binned
function [speed_array] = binWeibDist(weibVars, nMaxWindSpeed, nNumWeibBins)
    %-- From a tuple of <w_s(1)> <w_s(2)> which are <Lambda> and <k>
    %-- coefficients at a specific wind direction, binWeibDist() returns a
    %-- matrix speed_array where <s_a(1)> is the speed, <s_a(2)> is the
    %-- frequency that speed occurs.
    
    % Needed constants
    nBinDivisions = linspace(0, nMaxWindSpeed, nNumWeibBins+1);
    speed_array = zeros(nNumWeibBins,2);
    
    f = @(x) wblpdf(x,weibVars(1),weibVars(2));
    % Integrate and check
    %L = integral(f, 0,nMaxWindSpeed); % Should = 1
    for i = 1:nNumWeibBins
        speed_array(i,1) = (nBinDivisions(i) + nBinDivisions(i+1))/2;
        speed_array(i,2) = integral(f, nBinDivisions(i),nBinDivisions(i+1));
    end
    % Normalize
    speed_array(:,2) = speed_array(:,2)/sum(speed_array(:,2));
end