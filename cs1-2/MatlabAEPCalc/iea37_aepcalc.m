% Nick Baker
% Created for IEA Task 37 Wind Farm Layout Optimization Case Study 1
% Create 17 Sept 18

clear, close all;
% Needed for .yaml reading ability
addpath(genpath('/Users/nbaker/Documents/MATLAB/YAMLMatlab'));

rng('shuffle');        % Ensures we're using "random" numbers
%rng('default');         % For repeatable "random" numbers
nNumCycles = 2;

farmSize = 1;   % 0 = 9 turbines, 1 = 16 turbs, 2 = 36 turbs, 3 = 64 turbs
switch(farmSize)
    case 0
        fname_turb_loc = 'iea37-ex9.yaml';
        farm_rad = 900;
        plot_dimen = 1200;
        nNumRtrs = 9;
    case 1
        fname_turb_loc = 'iea37-ex16.yaml';
        farm_rad = 1300;
        plot_dimen = 1600;
        nNumRtrs = 16;
    case 2
        fname_turb_loc = 'iea37-ex36.yaml';
        farm_rad = 2000;
        plot_dimen = 2500;
        nNumRtrs = 36;
    case 3
        fname_turb_loc = 'iea37-ex64.yaml';
        farm_rad = 3000;
        plot_dimen = 3500;
        nNumRtrs = 64;
    otherwise
        error('Variable "FarmSize" not initilized properly');
end

% Get turbine data from .yaml
[~, fname_turb, fname_wr] = getTurbLocYAML(fname_turb_loc);
[turb_ci, turb_co, rated_ws, rated_pwr, turb_diam] = getTurbineYaml(fname_turb);
[wind_dir, wind_freq, wind_speed] = getWindRoseYaml(fname_wr);

% Initialize placeholders
OptFarm = zeros(nNumCycles,(nNumRtrs*2));
AEPlist = zeros(nNumCycles,1);     % format: [AEP, # function calls]
fnCallList = zeros(nNumCycles,1);  % format: [AEP, # function calls]
BestAEP = [1,0];                   % To store the location of our best AEP. (index, value)

%%{
% Run the optimization and save our results as we go
for i = 1:nNumCycles
    i                    % Output which iteration we're on
    x0 = ((rand(1,(nNumRtrs*2))-0.5)*farm_rad); % randomize turb locations between -farm radius and +farm radius
    [xopt,fopt,~,output] = optimize_iea37_wflocs(x0, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, farm_rad);
    
    OptFarm(i,:) = xopt;
    AEPlist(i) = fopt;           % Save the AEP
    output.funcCount
    fnCallList(i) = output.funcCount;    % Save the number of function calls
end
%%}

for j = 1:nNumCycles
    if (BestAEP(2) > AEPlist(i))    % If our new AEP is better (Remember negative switches)
        BestAEP(2) = AEPlist(i);    % Save it
        BestAEP(1) = i;             % And the index of which run we're on
    end
end
% Calculate AEP for best layout
best_coords = makeCoordStruct(OptFarm(BestAEP(1),:));
binned_AEP = calcAEP(best_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
AEPlist = [AEPlist,fnCallList];

% Save all values and write best result to a .yaml
%{
switch(farmSize)
    case 0
        csvwrite('ans-turbloc-opt9.csv',OptFarm)
        csvwrite('ans-aepcalls-opt9.csv',AEPlist)
        writeTurbLocYAML('iea37-baker-opt9.yaml', makeCoordArray(best_coords), fname_turb, fname_wr, binned_AEP);
    case 1
        csvwrite('ans-turbloc-opt16.csv',OptFarm)
        csvwrite('ans-aepcalls-opt16.csv',AEPlist)
        writeTurbLocYAML('iea37-baker-opt16.yaml', makeCoordArray(best_coords), fname_turb, fname_wr, binned_AEP);
    case 2
        csvwrite('ans-turbloc-opt36.csv',OptFarm)
        csvwrite('ans-aepcalls-opt36.csv',AEPlist)
        writeTurbLocYAML('iea37-baker-opt36.yaml', makeCoordArray(best_coords), fname_turb, fname_wr, binned_AEP);
    case 3
        csvwrite('ans-turbloc-opt64.csv',OptFarm)
        csvwrite('ans-aepcalls-opt64.csv',AEPlist)
        writeTurbLocYAML('iea37-baker-opt64.yaml', makeCoordArray(best_coords), fname_turb, fname_wr, binned_AEP);
    otherwise
        error('Variable "FarmSize" not initilized properly');
end
%}
% Plots the optimized rotor locations
bestIndex = BestAEP(1)
bestAEP = BestAEP(2)
%OptFarm(BestAEP(1),:)

%color_num = 0;  % 0 = blue, 1 = red, 2 = yellow
%plotFarm(best_coords, turb_diam, farm_rad, plot_dimen, color_num)
%}

%{
%--- Sample Output ---%
AEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
% Print AEP for each binned direction, with 5 digits behind the decimal.
AEP
% Print AEP summed for all directions
sum(AEP)
color_num = 3;  % 0 = blue, 1 = red, 2 = yellow
plotFarm(turb_coords, turb_diam, farm_rad, plot_dimen, color_num)
%[xopt, fopt, ~, ~] = optimize_iea37_wflocs(x0, nNumRtrs, turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr)
%}
function [xopt, fopt, exitflag, output] = optimize_iea37_wflocs(x0, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr, farm_rad)

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
        
        %--- Design Variables ---%
        turb_coords = makeCoordStruct(x);       % Rip to the struct for calcAEP() to work
        
        %-- Objective Function --%
        % AEP function
        dirAEP = calcAEP(turb_coords, wind_freq, wind_speed, wind_dir, turb_diam, turb_ci, turb_co, rated_ws, rated_pwr);
        scaleAEP = dirAEP / 1e6;    % scale down from MWh
        f = -sum(scaleAEP);
        
        %- Inequality Constraints -%
        % Constrain pair separation to be greater than two diameters from each other
        tc = [turb_coords.x,turb_coords.y];     % Strip the struct into a matrix
        cTurbSpace = pdist(tc, 'euclidean')';   % Get the distance between each turbine
        c = (2*turb_diam) - cTurbSpace;         % Constrain that the turbines are less than 2 diams apart
        
        % Constrain turbine locaton to be within windfarm boundary
        cTemp = hypot(turb_coords.x,turb_coords.y) - farm_rad; 
        c = [c; cTemp];                         % Append constraints to list
        
        %- Equality Constraints (None) -%
        ceq = [];
    end

    % ------------Call fmincon------------
    options = optimoptions(@fmincon, 'MaxFunctionEvaluations', 30000); % 'display', 'iter-detailed'); %,'MaxFunctionEvaluations',inf,'MaxIterations',Inf);
    [xopt, fopt, exitflag, output] = fmincon(@obj, x0, A, b, Aeq, beq, lb, ub, @con, options);
    
    
    % ------------Separate obj/con (do not change)------------
    function [f] = obj(x)
            [f, ~, ~] = objcon(x);
    end
    function [c, ceq] = con(x)
            [~, c, ceq] = objcon(x);
    end
end