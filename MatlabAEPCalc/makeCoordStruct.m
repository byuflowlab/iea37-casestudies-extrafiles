function [turb_coords] = makeCoordStruct(x0)
    % Takes coordinates in form turb_coords.x and .y and puts them into a
    % single array. For fmincon use.
    nNumRtrs = length(x0)/2;
    
    turb_coords.x = x0(1:nNumRtrs);
    turb_coords.y = x0((nNumRtrs+1):length(x0));
    %for i = 1:nNumRtrs
    %    turb_coords.x(i,1) = x0((i*2)-1); % x-coord
    %    turb_coords.y(i,1) = x0(i*2);     % y-coord
    %end
end