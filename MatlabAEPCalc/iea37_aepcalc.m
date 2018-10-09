% Nick Baker
% Created for IEA Task 37 Wind Farm Layout Optimization Case Study 1
% Create 17 Sept 18

clear all, close all
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

% Get turbine data from .yaml
[turb_coords, fname_turb, fname_wr] = getTurbLocYAML('iea37-ex16.yaml');
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);
%nNumRtrs = length(turb_coords.x);   % Pulls the number of turbines by how many x-coordinates we have.

%--- Sample Output ---%
AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
% Print AEP for each binned direction, with 5 digits behind the decimal.
AEP
% Print AEP summed for all directions
sum(AEP)
%[xopt, fopt, ~, ~] = optimize_iea37_wflocs(x0, nNumRtrs, turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)

function [xopt, fopt, exitflag, output] = optimize_iea37_wflocs(x0, nNumRtrs, turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)

    % ------------Starting point and bounds------------
    ub = [];
    lb = [];

    % ------------Linear constraints------------
    A = [];
    b = [];
    Aeq = [];
    beq = [];

    % ------------Objective and Non-linear Constraints------------
    function [f, c, ceq] = objcon(x)
        %--- Analysis Variables ---%
        [nNumPairs,~] = size(combnk(1:nNumRtrs,2));   % Calculates number of unique pairs for spacing constraints using combinatorics 
        
        %-- Objective Function --%
        % AEP function
        f = -calcAEP(x0, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
        
        %- Inequality Constraints -%
        % Constrain pair separation to be greater than two diameters from each other
        tc = [turb_coords.x,turb_coords.y];     % Strip the struct into a matrix
        c = (-(pdist(tc, 'euclidean')).^2 + ((2*turb_diam)^2))';
        
        % Constrain turbine locaton to be within windfarm boundary
        cTemp = -(hypot(turb_coords.x,turb_coords.y).^2) + (farm_rad)^2; 
        c = [c; cTemp]; % Append constraints to list
        
        %- Equality Constraints (None) -%
        ceq = [];
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'display', 'iter-detailed');
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end
end