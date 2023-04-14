function plotbar(data,overlap,color,xpos,width,transparency)
% add transparency (version200704)

if nargin < 6; transparency = 1; end
if nargin < 5; width = 0.7; end
if nargin < 4; xpos = 1:size(data,2); end
if nargin < 3; color = repmat([1,1,1],size(data,2),1); end
if nargin<2 % default: the same individual data will horizontally non-overlaply distributed
    overlap=1;
end

%% e.g.
% plotbar(data,1);
% ylim([0 100]);
% ax=gca;
% ax.XTick = 1:1:size(data,2);
% ax.XTickLabel = condlist(:);
% set(gca, 'YTick', 0:20:100);
% ylabel('Hit rate [%]');

% title('Gap detection performance');
% set(gca,'FontSize',12);
% % set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2 3]);
%
% yl = ylim;
% for x3 = 1:size(data,2)
%     p = mean(data(:,x3))*100;
%     y3 = yl(2)*0.1;%-0.1*(yl(2)-yl(1));
%     txt3 = [sprintf('%0.2f',p)];
%     text(x3,y3,txt3,'FontSize',12,'HorizontalAlignment','center');
% end
%
% filename = ['./Plots/' name '_behav_1hit'];
% print(filename,'-dpng','-r0');



hold on;

%% Mean bar plot
for i = 1:size(data,2)
    b1=bar(xpos(i), nanmean(data(:,i)), width);
    try
        b1.FaceColor=color(i,:); % transparent bar
    catch
        keyboard;
    end
    
    try
        b1.FaceAlpha = transparency;
    catch
        keyboard;
    end
    
    % b1.FaceColor=[1 1 0]; % yellow bar
    b1.EdgeColor= 'none';
    %     b1.LineWidth=2;
end
ax1 = gca;
hold(ax1, 'all');

x=[]; y=[];
for i=1:size(data,2)
    x=[x; xpos(i)*ones(size(data(:,i)))];
    y=[y; data(:,i)];
end

%% Individual subject data
% s1=scatter([ones(size(data(:,1)));2*ones(size(data(:,2)));3*ones(size(data(:,3)));4*ones(size(data(:,4)))]+.2,...
%     [data(:,1);data(:,2);data(:,3);data(:,4)],'filled', 'Parent', ax1);


if overlap
    for i=1:numel(x)
        idx=find(x==x(i) & y==y(i));
        if ~isempty(idx)
            for ii=1:numel(idx)
                if mod(ii,2)==0
                    %                     x(idx(ii))=x(idx(ii))+0.025*(ii-1);
                    x(idx(ii))=x(idx(ii))+0.05*(ii)/2;
                else
                    x(idx(ii))=x(idx(ii))-0.05*(ii-1)/2;
                end
            end
        end
    end
else
    x=x+.2;
end

% s1=scatter(x+.2,y,'filled', 'Parent', ax1);
s1=scatter(x,y,15,'o', 'Parent', ax1);
s1.MarkerEdgeColor = [.3 .3 .3];
s1.MarkerFaceColor = 'none';
% s1.MarkerFaceColor=[.5 .5 .5];
s1.LineWidth = 1;

%% Error bar
e1=errorbar(xpos,nanmean(data),nanstd(data)/sqrt(size(data,1)),'.');% SEM
% e1=errorbar(mean(data),std(data),'.');% Standard deviation
e1.Color=[0 0 0];
e1.LineWidth=2;
