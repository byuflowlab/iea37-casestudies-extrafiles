function testPlotLine(m,b)
    x = 9500:10:19000;
    y = m.* x + b;

    plot(x,y,'o')
    % Zero out and scale our upper limits
    ylim([0 2.25e4])
    xlim([0 1.9e4])
    % Make sure it's proportional
    yticks(linspace(0,2.25e4,5));
    xticks(linspace(0,1.9e4,5));
    gca.YRuler.Exponent = 0;
    axis square
end