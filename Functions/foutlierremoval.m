function y = foutlierremoval(x,numstd,nanmode)
if nargin<2; numstd=2; end
meanRT = nanmean(x);
stdRT = nanstd(x);

switch nanmode
    case 'nan'
        x(x>meanRT+numstd*stdRT) = nan; x(x<meanRT-numstd*stdRT) = nan;
    case 'remove'
        x(x>meanRT+numstd*stdRT) =[]; x(x<meanRT-numstd*stdRT) =[];
end

y=x;
end