% Combine the two tables for immediate and delayed alf

clear; close all; addpath(genpath('Functions'));

currentFolder = pwd;
[~,currentFolderName,~] = fileparts(pwd);
path_root = strrep(pwd,currentFolderName,''); % get the path of the parent folder

a1 = readtable(fullfile(path_root,'Project02_GiantProlific','alf_immediate_prolific','result_alf_immediate.csv'));
a2 = readtable(fullfile(path_root,'Project02_GiantProlific','alf_delayed_prolific','result_alf_delayed.csv'));

a = innerjoin(a1,a2,'Keys','participantID');

a.immediate_motivation_change = a.motivation_2_a1-a.motivation_1_a1; % 20230221 update fatigue rating
a.immediate_fatigue_change = a.fatigue_2_a1-a.fatigue_1_a1; % 20230221 update fatigue rating

a.motivation_change = a.motivation_1_a2-a.motivation_1_a1; % 20230221 update fatigue rating
a.fatigue_change = a.fatigue_1_a2-a.fatigue_1_a1; % 20230221 update fatigue rating

a.interval = seconds(a.starttime_a2 - a.endtime_a1)/60;


figure; clf;
yname = {'Number of correct words remembered after one presentation'};
yp = 'nCorrect_freerecall_1';
b = a.(yp);
% h = histogram(b,'Normalization','probability');
h = histfit(b,10,'kernel');
h(1).FaceColor= [168 115 211]/225;
h(1).EdgeColor = 'none';
h(2).Color = '#002147';
title(sprintf('Healthy Control (n = %d), median = %d',height(b),median(b)));
xlabel(yname); ylabel('Probability');
xlim([0,12]);
axis square;
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
exportgraphics(gcf,fullfile('Logbook',['hist_' yp '.png']),'Resolution',600);

figure; clf;
yname = {'Number of correct words remembered after learning'};
yp = 'nCorrect_freerecall_2';
b = a.(yp);
% h = histogram(b,'Normalization','probability');
h = histfit(b);
h(1).FaceColor= [168 115 211]/225;
h(1).EdgeColor = 'none';
h(2).Color = '#002147';
title(sprintf('Healthy Control (n = %d), median = %d',height(b),median(b)));
xlabel(yname); ylabel('Probability');
xlim([0,12]);
axis square;
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
exportgraphics(gcf,fullfile('Logbook',['hist_' yp '.png']),'Resolution',600);


figure; clf;
yname = {'Number of correct words recalled after 30 minutes'};
yp = 'nCorrect_freerecall_delayed';
b = a.(yp);
% h = histogram(b,'Normalization','probability');
h = histfit(b,12,'kernel');
h(1).FaceColor= [168 115 211]/225;
h(1).EdgeColor = 'none';
h(2).Color = '#002147';
title(sprintf('Healthy Control (n = %d), median = %d',height(b),median(b)));
xlabel(yname); ylabel('Probability');
xlim([0,12]);
axis square;
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
exportgraphics(gcf,fullfile('Logbook',['hist_' yp '.png']),'Resolution',600);




a.alf_forget = a.nCorrect_freerecall_delayed - a.nCorrect_freerecall_2;
figure; clf;
yname = {'Number of words forgotten in 30 minutes'};
yp = 'alf_forget';
b = a.(yp);
% h = histogram(b,'Normalization','probability');
h = histfit(b,10,'kernel');
h(1).FaceColor= [168 115 211]/225;
h(1).EdgeColor = 'none';
h(2).Color = '#002147';

hold on;
yl = ylim;
line([0 0],yl,'LineStyle','--','Color','k');
ylim(yl);
title(sprintf('Healthy Control (n = %d), median = %d',height(b),median(b)));
xlabel(yname); ylabel('Probability');
axis square;
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 5]);
exportgraphics(gcf,fullfile('Logbook',['hist_' yp '.png']),'Resolution',600);
% m = median(b)
% m = mean(b)

% duration of alf immediate
% figure; clf;
x = seconds(a1.endtime-a1.starttime)/60;
subplot(2,1,1);
histogram(x);
title('Length of ALF immediate');

x = seconds(a2.endtime-a2.starttime)/60;
subplot(2,1,2);
histogram(x);
title('Length of ALF delayed');
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);

a.ALF_forget_freerecall = a.nCorrect_freerecall_delayed - a.nCorrect_freerecall_2;
a.ALF_propremember_freerecall = a.nCorrect_freerecall_delayed./a.nCorrect_freerecall_2;
a.ALF_forget_recognise = a.propcorrect_recognise_delayed - a.propcorrect_recognise_2;
a = movevars(a,'ALF_forget_recognise','After','participantID');
a = movevars(a,'ALF_propremember_freerecall','After','participantID');
a = movevars(a,'ALF_forget_freerecall','After','participantID');
writetable(a,'result_alf.csv');



d = readtable(fullfile(path_root,'Project02_GiantProlific','demographics_oxcog_prolific','qualtrics','longCOVID_questionnaire.xlsx'));
a = innerjoin(a,d,'Keys','participantID');

d = readtable(fullfile(path_root,'Project02_GiantProlific','demographics_oxcog_prolific','demographic_prolific.csv'));
a = innerjoin(a,d,'Keys','participantID');

v = 'COVID_infection'; a.(v)(strcmp(a.(v),'Maybe')) = {'No'}; a = movevars(a,v,'After','participantID');
% v = 'COVID_infection'; a.(v)(strcmp(a.(v),'Maybe')) = {'Yes'}; a = movevars(a,v,'After','participantID');
v = 'brainfog'; a.(v)(strcmp(a.(v),'Maybe')) = {'Yes'}; a = movevars(a,v,'After','participantID');
v = 'memory_issue'; a.(v)(strcmp(a.(v),'Maybe')) = {'Yes'}; a = movevars(a,v,'After','participantID');
v = 'concentration_issue'; a.(v)(strcmp(a.(v),'Maybe')) = {'Yes'}; a = movevars(a,v,'After','participantID');




edges = [0,6,12,24,36];
groups = {'no infection','0-6m','6-12m','12-24m','24-36m'};
y = discretize(a.COVID_months,edges);
y = y + 1;
y(isnan(y)) = 1;
a.covidTimeGroup = groups(y)';
a = movevars(a,'covidTimeGroup','After','participantID');

% adjust all data with age
v = 'forget_freerecall';

writetable(a,'all_withALF.csv');


ois2 = readtable(fullfile('Results','result_ois2.csv'));
a = innerjoin(a,ois2,'Keys','participantID');
a.forget_ois2 = 1- a.OIS_DelayedObjectAccuracy./a.OIS_ImmediateObjectAccuracy;
a = movevars(a,'forget_ois2','After','participantID');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;clf;
param = 'memory_issue'; condlist = {'Yes','Maybe','No'};
% param = 'brainfog'; condlist = {'Yes','Maybe','No'};
xp = 'covidTimeGroup'; groups = {'no infection','0-6m','6-12m','12-24m','24-36m'};
% param = 'brainfog'; condlist = {'Yes','Maybe','No'};

red = [255 11 11];
redlight = [255 204 203];
blue = [38 17 249];

cm = [red;redlight;blue]/255;

ylist = {...
    'nCorrect_freerecall_1','nCorrect_freerecall_2','nCorrect_freerecall_delayed','ALF_propremember_freerecall',...
    'propcorrect_recognise_1','propcorrect_recognise_2','propcorrect_recognise_delayed','ALF_forget_recognise'...
    'immediate_motivation_change','immediate_fatigue_change','motivation_change','fatigue_change',};

for i = 1:length(ylist)
    yName = ylist{i};

    subplot(3,4,i);
    %     figure; clf;
    hold on;
    a.(xp) = categorical(a.(xp));

    a = movevars(a,param,'After','participantID');

    %     cm = cbrewer('qual','Set1',length(condlist));

    x = 1:length(groups);
    for k = 1:length(condlist)
        y = []; z = [];
        for s = x
            idx = find(a.(xp) == groups{s} & strcmp(a.(param),condlist{k}));
            data = a.(yName)(idx);
            data(isnan(data)) = [];
            y(s) = mean(data);
            z(s) = std(data)/sqrt(length(data));
            n(s) = length(data);

            txt = text(s+0.2,y(s),num2str(n(s)),'Color',cm(k,:),'FontWeight','bold','HorizontalAlignment','center');
        end

        eb = errorbar(x,y,z);
        eb.Color = cm(k,:);
        eb.LineWidth = 1;
        xlim([0.5 max(x)+0.5]);
        xticks(x);
        xticklabels(groups);

        grid on;
    end

    if i == 1
        legend(condlist,'Location','best');
    end
    %     title(strrep(param,'_',' '));
    ylabel(strrep(yName,'_',' '));

end
suptitle(strrep(param,'_',' '));
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);





xv = '0-6m';
yp = 'nCorrect_freerecall_2';
yes = a.(yp)(find(a.(xp) == xv & strcmp(a.(param),'Yes'))); mean(yes)
no =  a.(yp)(find(a.(xp) == xv & strcmp(a.(param),'No'))); mean(no)
[~,p] = ttest2(yes,no);
p


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;clf;
param = 'memory_issue'; condlist = {'Yes','Maybe','No'};
% param = 'brainfog'; condlist = {'Yes','Maybe','No'};
xp = 'covidTimeGroup'; groups = {'no infection','0-6m','6-12m','12-24m','24-36m'};
% param = 'brainfog'; condlist = {'Yes','Maybe','No'};

red = [255 11 11];
redlight = [255 204 203];
blue = [38 17 249];

cm = [red;redlight;blue]/255;

ylist = {'forget_ois2','DelayedSemanticAccuracy','ProportionCorrect_3item_4secs','Misbinding_Simple_3item_4secs',...
    'TMT_A','TMT_B','DSST_nCorrectResponse','AbsoluteError_1item_4secs'};

for i = 1:length(ylist)
    yName = ylist{i};

    subplot(3,4,i);
    %     figure; clf;
    hold on;
    a.(xp) = categorical(a.(xp));

    a = movevars(a,param,'After','participantID');

    %     cm = cbrewer('qual','Set1',length(condlist));

    x = 1:length(groups);
    for k = 1:length(condlist)
        y = []; z = [];
        for s = x
            idx = find(a.(xp) == groups{s} & strcmp(a.(param),condlist{k}));
            data = a.(yName)(idx);
            data(isnan(data)) = [];
            y(s) = mean(data);
            z(s) = std(data)/sqrt(length(data));
            n(s) = length(data);

            txt = text(s+0.2,y(s),num2str(n(s)),'Color',cm(k,:),'FontWeight','bold','HorizontalAlignment','center');
        end

        eb = errorbar(x,y,z);
        eb.Color = cm(k,:);
        eb.LineWidth = 1;
        xlim([0.5 max(x)+0.5]);
        xticks(x);
        xticklabels(groups);

        grid on;
    end

    if i == 1
        legend(condlist,'Location','best');
    end
    %     title(strrep(param,'_',' '));
    ylabel(strrep(yName,'_',' '));

end
suptitle(strrep(param,'_',' '));
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);




% age bin
figure;
hist(a.age);
edges = [20,30,40,50,60,80];
ageGroupName = {'20s','30s','40s','50s','>60'};
y = discretize(a.age,edges);
a.ageGroup = ageGroupName(y)';
a = movevars(a,'ageGroup','After','participantID');


% plot
param = 'memory_issue'; condlist = {'Yes','Maybe','No'};
% param = 'brainfog'; condlist = {'Yes','Maybe','No'};

red = [255 11 11];
redlight = [255 204 203];
blue = [38 17 249];

cm = [red;redlight;blue]/255;

figure;clf;
% ylist = {'forget_freerecall','forget_recognise','forget_ois2','Target_Simple_3item_4secs','DSST_nCorrectResponse','TMT_A','TMT_B','nCorrect_freerecall_1','n_repeated_freerecall','nCorrect_freerecall_delayed'};
% ylist = {'immediate_motivation_change','immediate_fatigue_change','nCorrect_freerecall_1','nCorrect_freerecall_delayed',...
%     'motivation_change','fatigue_change','forget_freerecall','forget_recognise'};
ylist = {...
    'motivation_1_a1','motivation_2_a1','motivation_1_a2','motivation_2_a2',...
    'fatigue_1_a1','fatigue_2_a1','fatigue_1_a2','fatigue_2_a2',...
    'confidence_1_a1','confidence_2_a1','confidence_2_a2','confidence_2_a2',...
    };
for i = 1:length(ylist)
    yName = ylist{i};

    subplot(3,4,i);
    %     figure; clf;
    hold on;
    a.ageGroup = categorical(a.ageGroup);

    a = movevars(a,param,'After','participantID');

    %     cm = cbrewer('qual','Set1',length(condlist));

    x = 1:length(ageGroupName);
    for k = 1:length(condlist)
        y = []; z = [];
        for s = x
            idx = find(a.ageGroup == ageGroupName{s} & strcmp(a.(param),condlist{k}));
            data = a.(yName)(idx);
            data(isnan(data)) = [];
            y(s) = mean(data);
            z(s) = std(data)/sqrt(length(data));
            n(s) = length(data);

            txt = text(s+0.2,y(s),num2str(n(s)),'Color',cm(k,:),'FontWeight','bold','HorizontalAlignment','center');
        end

        eb = errorbar(x,y,z);
        eb.Color = cm(k,:);
        eb.LineWidth = 1;
        xlim([0.5 max(x)+0.5]);
        xticks(x);
        xticklabels(ageGroupName);

        ylim([2,5]);
        yticks(2:0.5:5); grid on;
    end

    if i == 1
        legend(condlist,'Location','best');
    end
    %     title(strrep(param,'_',' '));
    ylabel(strrep(yName,'_',' '));

end
suptitle(strrep(param,'_',' '));
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);





figure;clf;
ylist = {...
    'nCorrect_freerecall_1','nCorrect_freerecall_2','nCorrect_freerecall_delayed','propremember_freerecall',...
    'propcorrect_recognise_1','propcorrect_recognise_2','propcorrect_recognise_delayed','forget_recognise'...
    'immediate_motivation_change','immediate_fatigue_change','motivation_change','fatigue_change',};

for i = 1:length(ylist)
    yName = ylist{i};

    subplot(3,4,i);
    %     figure; clf;
    hold on;
    a.ageGroup = categorical(a.ageGroup);

    a = movevars(a,param,'After','participantID');

    %     cm = cbrewer('qual','Set1',length(condlist));

    x = 1:length(ageGroupName);
    for k = 1:length(condlist)
        y = []; z = [];
        for s = x
            idx = find(a.ageGroup == ageGroupName{s} & strcmp(a.(param),condlist{k}));
            data = a.(yName)(idx);
            data(isnan(data)) = [];
            y(s) = mean(data);
            z(s) = std(data)/sqrt(length(data));
            n(s) = length(data);

            txt = text(s+0.2,y(s),num2str(n(s)),'Color',cm(k,:),'FontWeight','bold','HorizontalAlignment','center');
        end

        eb = errorbar(x,y,z);
        eb.Color = cm(k,:);
        eb.LineWidth = 1;
        xlim([0.5 max(x)+0.5]);
        xticks(x);
        xticklabels(ageGroupName);

        grid on;
    end

    if i == 1
        legend(condlist,'Location','best');
    end
    %     title(strrep(param,'_',' '));
    ylabel(strrep(yName,'_',' '));

end
suptitle(strrep(param,'_',' '));
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);



figure;clf;
ylist = {...
    'propcorrect_recognise_1','propcorrect_recognise_2','propcorrect_recognise_delayed'...
    'phit_recognise_1','phit_recognise_2','phit_recognise_delayed'...
    'pfa_novel_recognise_1','pfa_novel_recognise_2','pfa_novel_recognise_delayed'...
    'pfa_novel_recognise_1','pfa_old_recognise_2','pfa_old_recognise_delayed'...
    };

for i = 1:length(ylist)
    yName = ylist{i};

    subplot(4,3,i);
    %     figure; clf;
    hold on;
    a.ageGroup = categorical(a.ageGroup);

    a = movevars(a,param,'After','participantID');

    %     cm = cbrewer('qual','Set1',length(condlist));

    x = 1:length(ageGroupName);
    for k = 1:length(condlist)
        y = []; z = [];
        for s = x
            idx = find(a.ageGroup == ageGroupName{s} & strcmp(a.(param),condlist{k}));
            data = a.(yName)(idx);
            data(isnan(data)) = [];
            y(s) = mean(data);
            z(s) = std(data)/sqrt(length(data));
            n(s) = length(data);

            txt = text(s+0.2,y(s),num2str(n(s)),'Color',cm(k,:),'FontWeight','bold','HorizontalAlignment','center');
        end

        eb = errorbar(x,y,z);
        eb.Color = cm(k,:);
        eb.LineWidth = 1;
        xlim([0.5 max(x)+0.5]);
        xticks(x);
        xticklabels(ageGroupName);

        grid on;
    end

    if i == 1
        legend(condlist,'Location','best');
    end
    %     title(strrep(param,'_',' '));
    ylabel(strrep(yName,'_',' '));

end
suptitle(strrep(param,'_',' '));
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);


figure;clf;
% ylist = {'forget_freerecall','forget_recognise','forget_ois2','Target_Simple_3item_4secs','DSST_nCorrectResponse','TMT_A','TMT_B','nCorrect_freerecall_1','n_repeated_freerecall','nCorrect_freerecall_delayed'};
ylist = {'forget_ois2','DelayedSemanticAccuracy','ProportionCorrect_3item_4secs','Misbinding_Simple_3item_4secs',...
    'TMT_A','TMT_B','DSST_nCorrectResponse','AbsoluteError_1item_4secs'};
for i = 1:length(ylist)
    yName = ylist{i};

    subplot(2,4,i);
    %     figure; clf;
    hold on;
    a.ageGroup = categorical(a.ageGroup);

    a = movevars(a,param,'After','participantID');

    %     cm = cbrewer('qual','Set1',length(condlist));

    x = 1:length(ageGroupName);
    for k = 1:length(condlist)
        y = []; z = [];
        for s = x
            idx = find(a.ageGroup == ageGroupName{s} & strcmp(a.(param),condlist{k}));
            data = a.(yName)(idx);
            data(isnan(data)) = [];
            y(s) = mean(data);
            z(s) = std(data)/sqrt(length(data));
            n(s) = length(data);

            txt = text(s+0.2,y(s),num2str(n(s)),'Color',cm(k,:),'FontWeight','bold','HorizontalAlignment','center');
        end

        eb = errorbar(x,y,z);
        eb.Color = cm(k,:);
        eb.LineWidth = 1;
        xlim([0.5 max(x)+0.5]);
        xticks(x);
        xticklabels(ageGroupName);
    end

    if i == 1
        legend(condlist,'Location','best');
    end
    %     title(strrep(param,'_',' '));
    ylabel(strrep(yName,'_',' '));

end
suptitle(strrep(param,'_',' '));
set(findall(gcf,'-property','FontName'),'FontName','Diverda Sans Com');
% set(gcf,'PaperUnits','inches','PaperPosition',[0 0 5 10]);
