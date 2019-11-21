function [f,g,newDirsRad,newFreqs] = extrapolateFrequencies(oldFreqs,numNewDirs)
    %--- Given a list of frequencies around a windrose, splines and then reslices the trends.
    
    %-- Setup --%
    % Repeat the values on either side to get smooth slopes through endpoints
    oldFreqs3 = [oldFreqs,oldFreqs,oldFreqs];
    numOldDirs = length(oldFreqs);
    oldDirs = linspace((-2)*pi,4*pi-6*pi/(numOldDirs*3),(numOldDirs*3));
    
    %-- Make the spline --%
    fpp = spline(oldDirs,oldFreqs3);
    f = @(t) ppval(fpp,t);
    L = integral(f, 0,2*pi);
    
    %-- Normalize it for the right integral --%
    g = @(t) ppval(fpp,t)/L;
    %L2 = integral(g, 0,2*pi);   % L2 should = 1 (as a check)
    
    %-- Make new partitions --%
    %sizeBinRad = (2*pi)/numNewDirs;     % The radian size of each of our new bins
    %sizeHalfBin = sizeBinRad/2;         % Half that size, for math rasons
    %binBoundaries = linspace(-sizeHalfBin,(2*pi)-sizeHalfBin,(numNewDirs+1));
    %newDirsRad = linspace(0,2*pi-(2*pi/numNewDirs),numNewDirs);
    %newFreqs = zeros(numNewDirs,1);
    %for i=1:numNewDirs
    %    newFreqs(i) = integral(g, binBoundaries(i), binBoundaries(i+1));
    %end
    newDirsRad = linspace(0,2*pi-(2*pi/numNewDirs),numNewDirs);
    newFreqs = getFreqSlices(g, numNewDirs);
    
    %hold on
    %plot(oldDirs(13:24), oldFreqs)
    %plot(newDirs, f(newDirs));
    %plot(newDirs, newFreqs)
    %hold off
end