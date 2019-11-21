function [y] = calcWeibull(x, k, L)
    if (x >= 0)
        FirstPart = (k/L).*((x/L).^(k-1));
        ScndPart = exp(-((x/L).^k));
        y = FirstPart .* ScndPart;
    else
        y = 0;
    end 
end