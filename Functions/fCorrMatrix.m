function fCorrMatrix(A,param1,param2,corrType,pval,param1_name,param2_name)

if nargin<6
    param1_name = param1;
    param2_name = param2;
end

if nargin <5
    
    pval = 0.05;
end

% colourmap = cbrewer('seq','YlOrRd',10);
colourmap = cbrewer('seq','Oranges',10);
% colourmap_cold = cbrewer('seq','YlGnBu',10);
colourmap_cold = cbrewer('seq','Greens',10);
% figure(10000);clf; hold on;
% for i = 1:length(param1)
%     for j = 1:length(param2)
%         text(i,j,sprintf("%.2f",i/j));
%     end
% end
% xlim([0,length(param1)+1]);
% ylim([0,length(param1)+1]);

hold on;
for i = 1:length(param1)
    for j = 1:length(param2)
        a = A.(param1{i});
        b = A.(param2{j});
        
        I = find(isnan(a) | isnan(b));
        a(I) = [];
        b(I) = [];
        
        if ~isempty(a)
            [r,p] = corr(a, b, 'type',corrType);
            
            if r>0
                if r >= 0.9
                    c = colourmap(10,:);
                elseif r < 0.1
                    c = colourmap(1,:);
                elseif r < 0.2 && r>=0.1
                    c = colourmap(2,:);
                elseif r < 0.3 && r>=0.2
                    c = colourmap(3,:);
                elseif r < 0.4 && r>=0.3
                    c = colourmap(4,:);
                elseif r < 0.5 && r >=0.4
                    c = colourmap(5,:);
                elseif r < 0.6 && r>=0.5
                    c = colourmap(6,:);
                elseif r < 0.7 && r>=0.6
                    c = colourmap(7,:);
                elseif r < 0.8 && r>=0.7
                    c = colourmap(8,:);
                elseif r < 0.9 && r>=0.8
                    c = colourmap(9,:);
                end
            end
            if r<= 0
                if r <= -0.9
                    c = colourmap_cold(10,:);
                elseif r > -0.1
                    c = colourmap_cold(1,:);
                elseif r > -0.2 && r<=-0.1
                    c = colourmap_cold(2,:);
                elseif r > -0.3 && r<=-0.2
                    c = colourmap_cold(3,:);
                elseif r > -0.4 && r<=-0.3
                    c = colourmap_cold(4,:);
                elseif r > -0.5 && r <=-0.4
                    c = colourmap_cold(5,:);
                elseif r > -0.6 && r<=-0.5
                    c = colourmap_cold(6,:);
                elseif r > -0.7 && r<=-0.6
                    c = colourmap_cold(7,:);
                elseif r > -0.8 && r<=-0.7
                    c = colourmap_cold(8,:);
                elseif r > -0.9 && r<=-0.8
                    c = colourmap_cold(9,:);
                end
            end
            
            if isequal(param1,param2) % if param1 == param2
                if i/j<1
                    ifShow = 1;
                else
                    ifShow = 0;
                end
            else
                if isequal(param1{i},param2{j})
                    ifShow = 0;
                else
                    ifShow = 1;
                end
            end
            if ifShow
                scatter(i,j,1000,c,'filled','square');
                if p < pval
                    text(i,j,sprintf('%.2f',r),'FontWeight','bold','HorizontalAlignment','center');
                else
                    text(i,j,sprintf('%.2f',r),'FontAngle','italic','HorizontalAlignment','center');
                end
            end
        end
    end
    xlim([1-0.5 length(param1)+0.5]);
    ylim([1-0.5 length(param2)+0.5]);
end
if isequal(param1,param2) % if param1 == param2
    xlim([1-0.5 length(param1)-0.5]);
    ylim([1.5 length(param2)+0.5]);
end
xticks(1:1:length(param1_name));
xticklabels(strrep(param1_name,'_',' '));
xtickangle(45);
set(gca, 'XAxisLocation', 'top');


yticks(1:1:length(param2_name));
yticklabels(strrep(param2_name,'_',' '));
% ytickangle(45);
% set(gca, 'YAxisLocation', 'right');
hold off;

if length(param1) == length(param2)
    axis square;
    if length(param1)>6
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 length(param1)*0.55 length(param2)*0.55]);
    else
        set(gcf,'PaperUnits','inches','PaperPosition',[0 0 length(param1)*0.6 length(param2)*0.6]);
    end
else
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 length(param1)*0.5 length(param2)*1]);
end
set(findall(gcf,'-property','FontName'),'FontName','Roboto');
