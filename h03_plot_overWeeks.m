clear; close all;
q = readtable(fullfile('results','ect_02a_octal_participantID_matched.csv'));
q = sortrows(q,"weekID","ascend");

mkdir('figures');

% task = 'dsst';
% task = 'alf';
% task = 'omt';
% task = 'ois';
% task = 'tmt';
task = 'other';

switch task
    case 'other'
        plist = {'srt_RT','srt_cvRT','RT_go',...
            'RT_nogo','pfa_nogo','cvRT_go',...
            'crt_RT','crt_cvRT','crt_pcorrect',...
            };
        plistName = strrep(plist,'_',' ');
        n_plot = length(plist);
        nsubplotx = 3;
        nsubploty = 3;

    case 'omt'
        plist = {'OMT_ProportionCorrect','OMT_AbsoluteError',...
            'OMT_Imprecision','OMT_IdentificationTime','OMT_LocalisationTime',...
            'OMT_TargetDetection','OMT_Misbinding','OMT_Guessing'};
        plistName = {'Identification Accuracy','Location Error','Imprecision in localisation','Identification Time','Localisation Time','Target Detection Rate','Misbinding Rate','Guessing Rate'};
        n_plot = length(plist);
        nsubplotx = 3;
        nsubploty = 3;

    case 'alf'
        plist = {'ALF_nImmediateRecall','ALF_pImmediateRecall',...
            'ALF_nDelayedRecall','ALF_pDelayedRecall'};
        plistName = {'number of words immediately recalled','correct recognition','number of words recalled','correct recognition (delayed)'};
        n_plot = length(plist);
        nsubplotx = 2;
        nsubploty = 2;

    case 'ois'
        plist = {'ois_ImmediateObjectAccuracy','ois_ImmediateSemanticAccuracy','ois_ImmediateLocationError',...
            'ois_DelayedObjectAccuracy','ois_DelayedSemanticAccuracy','ois_DelayedLocationError'};
        plistName = {'Object identification accuracy','Semantic identification accuracy','Location error',...
            'Object identification accuracy','Semantic identification accuracy','Location error'};
        n_plot = length(plist);
        nsubplotx = 2;
        nsubploty = 3;
    case 'corsi'
        plist = {'CORSI_mean'};
        plistName = {'Mean location error'};
        n_plot = length(plist);
        nsubplotx = 1;
        nsubploty = 1;
    case 'rocf'
        plist = {'REY_copy_score','REY_recall_score'};
        plistName = {'Copy score','Immediate recall score'};
        n_plot = length(plist);
        nsubplotx = 2;
        nsubploty = 1;
    case 'dsst'
        plist = {'DSST_nCorrectResponse','DSST_pHit'};
        plistName = {'Number of correct responses','Accuracy'};
        n_plot = length(plist);
        nsubplotx = 2;
        nsubploty = 1;
    case 'tmt'
        plist = {'TMT_A','TMT_B'};
        plistName = {'TMT A','TMT B'};
        n_plot = length(plist);
        nsubplotx = 2;
        nsubploty = 1;
end



weeklist = 0:1:6;


% vname = 'DSST_nCorrectResponse';
for ii = 1:length(plist)
    subplot(nsubplotx,nsubploty,ii); hold on;
    variable = plist{ii};
    variableN = plistName{ii};

    x = nan(length(weeklist),1);
    y = nan(length(weeklist),1);
    z = nan(length(weeklist),1);
    for i = 1:length(weeklist)

        x(i) = weeklist(i);
        data = q.(variable)(q.weekID == weeklist(i)); data(isnan(data)) = [];
        y(i) = nanmean(data);
        z(i) = nanstd(data)/sqrt(length(data));

    end

    % figure(1);
    h = errorbar(x,y,z);
    h.Color = [0 0 0];
    h.LineWidth = 1;
    h.Marker = 'o';

    % plot(x,y,'*-');

    xticks(x);
    xlim([-0.5, max(x)+0.5]);
    xlabel('Week Number');
    ylabel(strrep(variable,'_',' '));
    title(variableN);
end

set(findall(gcf,'-property','FontName'),'FontName','Corbel','FontSize',9);
set(gcf,'PaperUnits','inches','PaperPosition',[0 0 2.8*nsubploty 2.3*nsubplotx]);
saveas(gcf,fullfile('figures',['figure1_byWeek_' task '.png']));
