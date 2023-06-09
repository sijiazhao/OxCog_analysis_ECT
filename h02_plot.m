clear; close all;

q = readtable('sart_vs_garden.csv');

mkdir('figures');


%% PLOT BAR
data = q;
variable = 'medianRT';
variableName = 'mediant RT (s)';
h02_plot_extra;

data = q;
variable = 'cvRT';
variableName = 'coefficient of variation';
h02_plot_extra;

data = q;
variable = 'pfa';
variableName = 'error of commission rate';
h02_plot_extra;











data = q;
variable = 'fatigue_accumulation';
variableName = 'accumulated fatigue';
condlist = {'sart','garden'};
plist = {[variable '_sart'], [variable '_garden']};
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

% [~,p] = ttest(a(:,2),a(:,3));
% if p<0.001; str = '***';
% elseif p<0.01; str = '**';
% elseif p<0.05; str = '*';
% else; str = 'n.s.';end
% text(2.5,yl(1)+0.1*mean(yl),str,'FontSize',15,'HorizontalAlignment','center','VerticalAlignment','bottom');

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
saveas(gcf,fullfile('figures',['03_bar_' variable '.png']));

