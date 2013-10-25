function [fl,fu] = calcCutOff(cf,bw);

if isnan(bw)
    fl=NaN;
    fu=NaN;
else
    fl = floor(cf / 2.^(bw./2));
    fu = ceil(cf * 2.^(bw./2));
end

return;