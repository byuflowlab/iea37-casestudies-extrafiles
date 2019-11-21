% Nick Baker
% 03 Apr 19

function [] = plotWindRoseFreq(binnedDirsRad, binnedFreqs)
    % To plot the Wind Frequency distribution
    binnedDirsRad =  [binnedDirsRad,binnedDirsRad(1)];
    binnedFreqs = [binnedFreqs,binnedFreqs(1)];
    polarplot(binnedDirsRad,binnedFreqs,'LineWidth',2);   % Plot outside edge
    ax = gca;
    ax.FontSize =8;
    ax.ThetaZeroLocation = 'top';
    ax.ThetaDir = 'clockwise';
    ax.ThetaAxisUnits = 'degrees';
    rticks([.0025 .0075 .0125])
    rticklabels({'.25%','.75%','1.25%'})
end