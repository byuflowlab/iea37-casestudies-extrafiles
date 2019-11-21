function [newFreqs] = getFreqSlices(g, numNewDirs)
    
    %-- Make new partitions --%
    sizeBinRad = (2*pi)/numNewDirs;     % The radian size of each of our new bins
    sizeHalfBin = sizeBinRad/2;         % Half that size, for math rasons
    binBoundaries = linspace(-sizeHalfBin,(2*pi)-sizeHalfBin,(numNewDirs+1));
    newFreqs = zeros(numNewDirs,1);
    for i=1:numNewDirs
        newFreqs(i) = integral(g, binBoundaries(i), binBoundaries(i+1));
    end
end