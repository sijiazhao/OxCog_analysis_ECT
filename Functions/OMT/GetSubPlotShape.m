function [r, c] = GetSubPlotShape(n)
% function [r, c] = GetSubPlotShape(n)
% use ceil&sqrt to get good shape for subplots depending on number
% will be rectangular, roughly square.
% If only one argument out (e.g. [r] = GetSubplotShape(n);) then it will be
% a 1x2 vector of [r c]

r = ceil(sqrt(n));
c = ceil(n / r);
end