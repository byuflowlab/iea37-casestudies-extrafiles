function [] = plotSubfigLgd(lgdColor,lgdName,nNumPar)
    addpath(genpath('/Users/nbaker/Documents/MATLAB/subplot_tight/subplot_tight'));
    
    if (nNumPar == 5)
        [m,n,~] = getSubfigPos((nNumPar+1), nNumPar);   % Get the legend positioning
        p = 2;
    else
        [m,n,p] = getSubfigPos((nNumPar+1), nNumPar);   % Get the legend positioning
    end
    
    hL = subplot_tight(m,n,p);                            % Plot the legend
    poshL = get(hL,'position');     % Getting its position
    lgd = legend(hL,lgdColor,lgdName);
    set(lgd,'position',poshL);                      % Adjusting legend's position
    %axis(hL,'off');                                 % Turning its axis off
    grid off;
    set(gca,'xtick',[])
    set(gca,'ytick',[])
    set(gca,'XColor','none')
    set(gca,'YColor','none')
    pbaspect([1 1 1])                     % Make the plots actually square
end