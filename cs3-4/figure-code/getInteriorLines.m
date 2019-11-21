function [newBndryPts] = getInteriorLines(oldBndryEqs, oldBndryPts, offset)
    %-- Given a set of lines forming a closed polygon, adjusts each line by
    % <offset> towards the interior
    nNumLines = length(oldBndryEqs);
    bUpInside = true(nNumLines,1);       % 0 = down, 1 = up
    VectorAngle = zeros(nNumLines,1);
    newBndryEqs = oldBndryEqs;
    newBndryPts = zeros(length(oldBndryPts),2);
    %figure(1)
    %plotBorsseleBoundary(oldBndryPts)  % Debug check
    
    %-- Find the interior --%
    %- Loop clockwise through points -%
    for i = 1:(nNumLines-1)
        %slope(i) = (oldBndryPts(i+1,2)-oldBndryPts(i,2))/(oldBndryPts(i+1,1)-oldBndryPts(i,1));
        
        VectorAngle(i) = atan2(oldBndryPts(i+1,2)-oldBndryPts(i,2),oldBndryPts(i+1,1)-oldBndryPts(i,1));
        NormalAngle = (VectorAngle(i) - (pi/2));   % rotate 90deg clockwise
        %- Make sure we're between -pi and pi -%
        if (NormalAngle < -pi)                  % If we passed the axis
            NormalAngle = NormalAngle + (2*pi); % Get us right
        end
        %- Determine if shift should be up or down -%
        if (NormalAngle <= 0)
            bUpInside(i) = false;
        end
    end
    % Do line from last point to first point
    VectorAngle(nNumLines) = atan2(oldBndryPts(1,2)-oldBndryPts(nNumLines,2),oldBndryPts(1,1)-oldBndryPts(nNumLines,1));
    NormalAngle = (VectorAngle(nNumLines) - (pi/2));   % rotate 90deg clockwise
    %- Make sure we're between -pi and pi -%
    if (NormalAngle < -pi)                  % If we passed the axis
        NormalAngle = NormalAngle + (2*pi); % Get us right
    end
    %- Determine if shift should be up or down -%
    if (NormalAngle <= 0)
        bUpInside(nNumLines) = false;
    end
    
    %-- Move the correct direction (up or down) --%s
    for i = 1:nNumLines
        if (bUpInside(i))
            newBndryEqs(i,2) = oldBndryEqs(i,2) + offset*sqrt(1+oldBndryEqs(i,1)^2);
        else
            newBndryEqs(i,2) = oldBndryEqs(i,2) - offset*sqrt(1+oldBndryEqs(i,1)^2);
        end
    end
    
    %-- Find intersection points--%
    for i = 1:(nNumLines-1)
        % x = (b1-b2)/(m2-m1)
        newBndryPts(i+1,1) = (newBndryEqs(i,2) - newBndryEqs(i+1,2))/(newBndryEqs(i+1,1) - newBndryEqs(i,1));
        % y = m*x + b
        newBndryPts(i+1,2) = newBndryEqs(i,1)*newBndryPts(i+1,1) + newBndryEqs(i,2);
    end
    
    %- Do first point -%
    newBndryPts(1,1) = (newBndryEqs(nNumLines,2) - newBndryEqs(1,2))/(newBndryEqs(1,1) - newBndryEqs(nNumLines,1));        
    newBndryPts(1,2) = newBndryEqs(1,1)*newBndryPts(1,1) + newBndryEqs(1,2);
    %figure(2)
    %plotBorsseleBoundary(newBndryPts)  % Debug check
end