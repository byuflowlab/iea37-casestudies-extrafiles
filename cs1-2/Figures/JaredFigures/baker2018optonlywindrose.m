% Nicholas F. Baker
% Benchmark Competition 1
%    Jensen Cosine implementation
% Begun 7 Mar 2018
% Wind rose figure code

nNumRtrs = 9; % Number of rotors in field
nDegBucket = 1;      % How many degrees to have in each calculated bucket
nLastDeg = 360 - nDegBucket;
% Vars for Frequency calc
degDir = (0:nDegBucket:nLastDeg)';       % A list of every direction to calculate Frequency at
nNumBuckets = length(degDir);
FreqWindDir = zeros(nNumBuckets,1);
[nNumPairs,~] = size(combnk(1:nNumRtrs,2));   % Calculates number of unique pairs for spacing constraints using combinatorics 

nCntr = 0;
for i = 0:nDegBucket:nLastDeg
    nCntr = nCntr + 1;
    FreqWindDir(nCntr) = BenchmarkWindDist(i);
end 
CheckRose = sum(FreqWindDir);

% To plot the Wind Frequency distribution
thetaPlot = deg2rad(0:nDegBucket:(360-nDegBucket));
polarplot(thetaPlot,FreqWindDir)
ax = gca;
ax.ThetaZeroLocation = 'top';
ax.ThetaDir = 'clockwise';

function [DistValue] = BenchmarkWindDist(degDir)
    % Given a wind diection (degDir), returns the frequency, as given in
    % Benchmark paper
    
    % Vars from Frequency calc
    Mu(1) = 180;       % (deg) Specified in benchmark paper.
    Mu(2) = 350;       % (deg) Specified in benchmark paper.
    Mu(3) = -10;       % (deg) Specified in benchmark paper.
    Sig(1) = 20;      % (deg) Specified in benchmark paper.
    Sig(2) = 40;      % (deg) Specified in benchmark paper.
    W(1) = 0.5;    % coefficient to weight frequency. Specified in benchmark paper.
    W(2) = 0.5;    % coefficient to weight frequency. Specified in benchmark paper.
    
    FreqWindDir1stTerm = W(1)*( sqrt(1/ (2*pi*(Sig(1)^2)) ) ) * exp( -( ((degDir-Mu(1)).^2)./(2.*Sig(1)^2) ) );
    FreqWindDir2ndTerm = W(2)*( sqrt(1/ (2*pi*(Sig(2)^2)) ) ) * exp( -( ((degDir-Mu(2)).^2)./(2.*Sig(2)^2) ) );
    FreqWindDir3rdTerm = W(2)*( sqrt(1/ (2*pi*(Sig(2)^2)) ) ) * exp( -( ((degDir-Mu(3)).^2)./(2.*Sig(2)^2) ) );
    DistValue = FreqWindDir1stTerm + FreqWindDir2ndTerm + FreqWindDir3rdTerm;
end