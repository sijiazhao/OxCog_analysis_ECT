function fPlotScatterCorr(x,y,z,corrType)

msize = 20;
% corrPlotMode = 'rank';
corrPlotMode = 'value';

% x(x==max(x)) = nan;
% x(x==max(x)) = nan;

if ~isempty(z)
    % x = nanzscore(x);
    M1 = fitlm(z,x);
    rX = M1.Residuals.Raw;
    
    % y = nanzscore(y);
    M2 = fitlm(z,y);
    rY = M2.Residuals.Raw;
    
    x = rX;
    y = rY;
end

temp = [x,y];
I = any(isnan(temp), 2);
x(I) = [];
y(I) = [];

switch corrPlotMode
    case 'rank'
        % -- Plot rank
        Data = x; [~,p] = sort(Data,'descend'); r = 1:length(Data); r(p) = r; a = r';
        Data = y; [~,p] = sort(Data,'descend'); r = 1:length(Data); r(p) = r; b = r';
    otherwise
        a = x;
        b = y;
end

s = scatter(a,b,msize,'o');
s.MarkerFaceColor = 'b';
% s.MarkerFaceColor = 'none';
s.MarkerEdgeColor = 'b';
s.MarkerFaceAlpha = 0.2;

p = polyfit(a,b,1);
f = polyval(p,a);
plot(a,f,'Color','k','LineWidth',2);
[r,p] = corr(a, b, 'type',corrType);
xl = xlim;
yl = ylim;

% p = p*nBonferroni;
if p < 0.001
    str = sprintf('%s r = %.2f, p < 0.001',corrType,r);
else
    str = sprintf('%s r = %.2f, p = %.3f',corrType,r,p);
end
title(str,'FontSize',9);


% xticks(0:10:100);
% % yticks(-90:25:65);
% yticks(-20:20:80);

xlim(xl);
ylim(yl);
axis square;

