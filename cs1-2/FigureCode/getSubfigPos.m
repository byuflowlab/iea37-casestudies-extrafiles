% Nick Baker
% 2 Dec 18
% If a figure of subfigs, finds where to put the current participant.

function [m,n, nParNum] = getSubfigPos(nParNum, nNumPlt)
    nNumPlt = nNumPlt + 1;      % Add a plot for the legend plot
    
    % Number of columns to have
    if nNumPlt > 6      % If we have more than 6 figures
        n = 3;            % Make three columns
    elseif nNumPlt > 2 
        n = 2;
    else
        n = 1;
    end
    
    % Number of rows to have
    if nNumPlt > 9      % If we have more than 6 figures
        m = 4;            % Make three columns
    elseif nNumPlt > 4
        m = 3;
    elseif nNumPlt > 1
        m = 2;
    else
        m = 1;
    end
    
    if (nNumPlt == 6) && (nParNum > 1)
        nParNum = nParNum + 1;
    end
end