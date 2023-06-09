% condlist = {'short sart','sart','garden','short garden'};
% plist = {[variable '_77_sart'], [variable '_sart'], [variable '_garden'], [variable '_77_garden']};
% cmap = [252, 140,090; 255,223,146;144,190,224;075,116,178]/255;

condlist = {'sart','garden','short garden'};
plist = {[variable '_sart'], [variable '_garden'], [variable '_77_garden']};
cmap = [255,223,146;144,190,224;075,116,178]/255;

a = [];
for k = 1:length(plist)
    a = [a, data.(plist{k})];
end


% bar plot
figure; clf;
hold on;

x = 1:length(condlist);
y = mean(a);
z = std(a)/sqrt(height(a));
for i = 1:length(x)
    h = bar(x(i),y(i));
    h.FaceColor = cmap(i,:);
    h.EdgeColor = 'none';

    str = sprintf('%.2f',y(i));
    text(x(i),0,str,'FontSize',15,'HorizontalAlignment','center','VerticalAlignment','bottom');
end

x1 = x; x1(1) = x1(1) + 0.25; x1(length(x1)) = x1(length(x1)) - 0.25;
for i = 1:height(a)
    plot(x1,a(i,:),'Color',[0.5 0.5 0.5]);
end

eb = errorbar(x,y,z);
eb.Color = 'k';
eb.LineStyle = 'none';
eb.LineWidth = 2;
eb.CapSize = 15;

xticks(x); xticklabels(condlist);
% ylim([0,50]);
axis square;

[~,p] = ttest(a(:,1),a(:,2));
if p<0.001; str = '***';
elseif p<0.01; str = '**';
elseif p<0.05; str = '*';
else; str = 'n.s.';end
yl = ylim;
text(1.5,yl(1)+0.1*mean(yl),str,'FontSize',15,'HorizontalAlignment','center','VerticalAlignment','bottom');

[~,p] = ttest(a(:,2),a(:,3));
if p<0.001; str = '***';
elseif p<0.01; str = '**';
elseif p<0.05; str = '*';
else; str = 'n.s.';end
text(2.5,yl(1)+0.1*mean(yl),str,'FontSize',15,'HorizontalAlignment','center','VerticalAlignment','bottom');

% [~,p] = ttest(a(:,3),a(:,4));
% if p<0.001; str = '***';
% elseif p<0.01; str = '**';
% elseif p<0.05; str = '*';
% else; str = 'n.s.';end
% text(3.5,yl(1)+0.1*mean(yl),str,'FontSize',15,'HorizontalAlignment','center','VerticalAlignment','bottom');

title(sprintf('n=%d',height(data)));
ylabel(variableName);
set(findall(gcf,'-property','FontName'),'FontName','Avenir LT Pro 45 Book');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
saveas(gcf,fullfile('figures',['01_bar_' variable '_77.png']));


%% PLOT CORRELATION
figure; hold on;
a = [];
plotcondition = [1,2];
% plotcondition = [2,3];
for k = plotcondition
    a = [a, data.(plist{k})];
end
x = a(:,1);
y = a(:,2);
defaultColor = [0.7 0.7 0.7];
type = 'Pearson';

sc = scatter(x,y);
sc.SizeData = 50;
sc.MarkerFaceColor = defaultColor;
sc.MarkerEdgeColor = 'none';
sc.MarkerFaceAlpha = 1;


ri = isnan(x)| isnan(y);
x(ri) = []; y(ri) = [];
p = polyfit(x,y,1);
f = polyval(p,x);
plot(x,f,'Color',[0 0 0],'LineWidth',1,'LineStyle','--');

[r,p] = corr(x,y,'type',type);
if p<0.001
    txt= sprintf('r=%.2f, p<0.001 (n=%d)',r,length(x));
else
    txt= sprintf('r=%.2f, p=%.3f (n=%d)',r,p,length(x));
end

axis square;

xlim([nanmean(x)-2*nanstd(x),nanmean(x)+2*nanstd(x)]);
ylim([nanmean(y)-2*nanstd(y),nanmean(y)+2*nanstd(y)]);

xlabel(condlist{plotcondition(1)});
ylabel(condlist{plotcondition(2)});

txt = {variableName; txt};
tt = title(txt);
tt.FontWeight = 'normal';
tt.FontSize = 10;
tt.Color = [0 0 0];
tt.HorizontalAlignment = "center";

set(findall(gcf,'-property','FontName'),'FontName','Avenir LT Pro 45 Book');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
saveas(gcf,fullfile('figures',['02_corr_' variable '.png']));
