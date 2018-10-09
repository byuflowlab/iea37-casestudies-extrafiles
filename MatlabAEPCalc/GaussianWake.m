function [loss] = GaussianWake(frame_coords, turb_diam)
    %-- Return each turbine's total loss due to wake from upstream turbines
    % Equations and values explained in <iea37-wakemodel.pdf>
    num_turb = length(frame_coords.x);
    
    % Array holding the wake deficit seen at each turbine
    loss = zeros(num_turb,1);

    for i = 1:num_turb                   % Looking at each turb (Primary)
        loss_array = zeros(num_turb,1);  % Calculate the loss from all others
        for j = 1:num_turb               % Looking at all other turbs (Target)
            x = frame_coords.x(i) - frame_coords.x(j);   % Calculate the x-dist
            y = frame_coords.y(i) - frame_coords.y(j);   % And the y-offset
            
            loss_array(j) = GWakeEq(x, y, turb_diam);
            % Total wake losses from all upstream turbs, using sqrt of sum of sqrs
            loss(i) = sqrt(sum(loss_array.^2));
        end
    end
end