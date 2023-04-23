clear;

q = readtable('all.csv');
a = readtable('result_alf.csv');
a.starttime_a1 = [];
a.endtime_a1 = [];
a.starttime_a2 = [];
a.endtime_a2 = [];
q = mergetable(q,a);

currentFolder = pwd;
[~,currentFolderName,~] = fileparts(pwd);
path_root = strrep(pwd,currentFolderName,''); % get the path of the parent folder


d = readtable(fullfile(path_root,'Project02_GiantProlific','demographics_oxcog_prolific','qualtrics','longCOVID_questionnaire.xlsx'));
q = mergetable(q,d);

d = readtable(fullfile(path_root,'Project02_GiantProlific','demographics_oxcog_prolific','demographic_prolific.csv'));
q = mergetable(q,d);

v = 'COVID_infection'; q.(v)(strcmp(q.(v),'Maybe')) = {'No'}; q = movevars(q,v,'After','participantID');
% v = 'COVID_infection'; a.(v)(strcmp(a.(v),'Maybe')) = {'Yes'}; a = movevars(a,v,'After','participantID');
v = 'brainfog'; q.(v)(strcmp(q.(v),'Maybe')) = {'Yes'}; q = movevars(q,v,'After','participantID');
v = 'memory_issue'; q.(v)(strcmp(q.(v),'Maybe')) = {'Yes'}; q = movevars(q,v,'After','participantID');
v = 'concentration_issue'; q.(v)(strcmp(q.(v),'Maybe')) = {'Yes'}; q = movevars(q,v,'After','participantID');



% x = table2array(q(:,{'OIS_ImmediateLocationError','OIS_DelayedLocationError','OMT_AbsoluteError','REY_recall_score','OIS_ImmediateObjectAccuracy','OIS_DelayedObjectAccuracy','OMT_ProportionCorrect',...
%     'ALF_propremember_freerecall','ALF_forget_recognise'}));

x = table2array(q(:,{'REY_recall_score','OIS_ImmediateObjectAccuracy','OIS_DelayedObjectAccuracy','OMT_ProportionCorrect',...
    'ALF_propremember_freerecall','ALF_forget_recognise'}));
% x = table2array(q(:,{'REY_recall_score','OIS_ImmediateObjectAccuracy','OIS_DelayedObjectAccuracy','OIS_ImmediateLocationError','OIS_DelayedLocationError','OMT_ProportionCorrect','OMT_AbsoluteError'}));
[~,score] = pca(x);
q.memory_pca = score(:,1);

x = table2array(q(:,{'OIS_ImmediateObjectAccuracy','OIS_ImmediateLocationError','OMT_ProportionCorrect','OMT_AbsoluteError','nCorrect_freerecall_1'}));
[~,score] = pca(x);
q.memory_pca1 = score(:,1);

x = table2array(q(:,{'REY_recall_score','OIS_ImmediateLocationError','OIS_DelayedLocationError','OMT_AbsoluteError',}));
[~,score] = pca(x);
q.memory_pca1b = score(:,1);

x = table2array(q(:,{'ALF_propremember_freerecall','ALF_forget_recognise','OIS_DelayedObjectAccuracy','OIS_DelayedLocationError'}));
[~,score] = pca(x);
q.memory_pca2 = score(:,1);

x = table2array(q(:,{'DSST_nCorrectResponse','TMT_A','TMT_B'}));
[~,score] = pca(x);
q.ef_pca = score(:,1);

writetable(q,'all_withALF_withLongCOVID.csv');


function q = mergetable(q,d)
colname = d.Properties.VariableNames;
for j = 2:length(colname)
    if iscell(d.(colname{j})(1))
        %         q.(colname{j}) = cell(height(q),1);
        q.(colname{j}) = repmat({''},height(q),1);
    else
        q.(colname{j}) = nan(height(q),1);
    end

    for i = 1:height(q)

        idx = find(strcmp(d.participantID,q.participantID{i}));
        if ~isempty(idx)
            q.(colname{j})(i) = d.(colname{j})(idx);
        end
    end
end
end
