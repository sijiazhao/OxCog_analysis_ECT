function y = sillyCell2num(x)
% STUPID PSYCHOJS BUG
% DONT UNDERSTAND WHY THERE IS INVISIBLE CHARACTER IN FRONT OF THE FIRST
% COLUMN, LET ME MANUALLY REMOVE IT
y = nan(size(x));
% numel(x)
for i = 1:numel(x)
    x1 = x{i};
    x1(1) = [];
    x1 = str2num(x1);
    if isempty(x1)
        return;
    else
        y(i) = x1;
    end
end