function [loss_perc] = GWakeEq(x, y, turb_diam)
    % Constant thrust coefficient
    CT = 4.0*1./3.*(1.0-1./3.);
    % Constant, relating to a turbulence intensity of 0.075
    k = 0.0324555;
    loss_perc = 0;      % Default to no loss
    
    if (x > 0)          % If Primary is downwind of the Target
        sigma = k*x + turb_diam/sqrt(8);         % Calculate the wake loss
        % Simplified Bastankhah Gaussian wake model
        exponent = -0.5 * (y/sigma)^2;
        radical = 1 - CT/(8*sigma^2 / turb_diam^2);
        loss_perc = (1.-sqrt(radical)) * exp(exponent);
    end
    % Note that if the Target is upstream, loss is defaulted to zero
end