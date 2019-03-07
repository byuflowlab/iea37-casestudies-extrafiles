function [x0] = makeCoordArray(turb_coords)
    % Takes coordinates in form turb_coords.x and .y and puts them into a
    % single array. For fmincon use.
    %nNumRtrs = length(turb_coords.x);
    %x0 = zeros(nNumRtrs*2,1);
    
    x0 = [turb_coords.x, turb_coords.y];
    %for i = 1:nNumRtrs
    %    x0((i*2)-1) = turb_coords.x(i); % x-coord
    %    x0(i*2) = turb_coords.y(i);     % y-coord
    %end
end