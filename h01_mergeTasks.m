close all; clear;
try
    addpath(genpath('E:\GitHub\Functions'));
end

if ismac
    path_root = '/Users/sijiazhao/Documents/GitHub/Masud_OxfordCognition/ECT/';
else
    path_root = 'E:\GitHub\Masud_OxfordCognition\Project13_ECT';
end

mkdir('results');

alf_immediate = readtable(fullfile(path_root, 'ect_alf', 'results', 'result_alf_immediate.csv')); alf_immediate = handleMultiEntry(alf_immediate);
alf_delayed = readtable(fullfile(path_root, 'ect_alf_delayed', 'results', 'result_alf_delayed.csv')); alf_delayed = handleMultiEntry(alf_delayed);
alf_delayed.ALF_nDelayedRecall = alf_delayed.nCorrect_freerecall_delayed;
alf_delayed.ALF_pDelayedRecall = alf_delayed.propcorrect_recognise_delayed;
alf_immediate.ALF_nImmediateRecall = alf_immediate.nCorrect_freerecall_1;
alf_immediate.ALF_pImmediateRecall = alf_immediate.propcorrect_recognise_1;

ois = readtable(fullfile(path_root, 'ect_ois', 'results', 'ois.csv')); ois = handleMultiEntry(ois);
omt = readtable(fullfile(path_root, 'ect_omt', 'results', 'omt.csv')); omt = handleMultiEntry(omt);

tmt = readtable(fullfile(path_root, 'ect_tmt', 'results', 'a.csv')); tmt = handleMultiEntry(tmt);
dsst = readtable(fullfile(path_root, 'ect_dsst', 'results', 'a.csv')); dsst = handleMultiEntry(dsst);

srt = readtable(fullfile(path_root, 'ect_srt', 'results', 'a.csv')); srt = handleMultiEntry(srt);
crt = readtable(fullfile(path_root, 'ect_crt', 'results', 'a.csv')); crt = handleMultiEntry(crt);
gonogo = readtable(fullfile(path_root, 'ect_gonogo', 'results', 'a.csv')); gonogo = handleMultiEntry(gonogo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge all tables
taskList = {'dsst','alf_immediate','alf_delayed','ois','omt','tmt','srt','crt','gonogo'};
% taskList = {'dsst','alf_immediate','alf_delayed','ois','omt','tmt'};

if 0 % IF ONLY KEEP SUBJECTS WHO COMPLETED EVERYTHING
    alf = innerjoin(alf_immediate,alf_delayed,'LeftKeys',1,'RightKeys',1);
    q = alf;
    % q = innerjoin(q,alf,'LeftKeys',1,'RightKeys',1);
    for i = 3:length(taskList)
        string = ['q = innerjoin(q,' taskList{i} ',''LeftKeys'',1,''RightKeys'',1);'];
        eval(string);
    end
else
    q = table;

    for j = 1:length(taskList)
        eval(['x = ' taskList{j} ';']);
        t = readtable(fullfile('ParameterNameList',['parameterList_task_' taskList{j} '.xlsx']));
        plist = t.plist;
        typelist = t.typelist;

        if j == 1
            q.entryID = x.entryID;
            q.participantID = x.participantID;
        end

        for i  = 1:length(plist)
            switch typelist{i}
                case 'num'
                    q.(plist{i}) = nan(size(q.entryID));
                case 'time'
                    q.(plist{i}) = NaT(size(q.entryID));
                case 'cell'
                    q.(plist{i}) = cell(size(q.entryID));
            end
        end

        for s = 1:size(x,1)
            idx = find(strcmp(q.entryID, x.entryID{s}));
            if ~isempty(idx)
                q.participantID(idx) = x.participantID(s);
                for i = 1:length(plist)

                    switch typelist{i}
                        case 'num'
                            q.(plist{i})(idx) = x.(plist{i})(s);
                        case 'time'
                            q.(plist{i})(idx) = x.(plist{i})(s);
                        case 'cell'
                            q.(plist{i}){idx} = x.(plist{i}){s};
                    end
                end

            else
                idx2 = size(q,1) + 1;
                columnName = q.Properties.VariableNames;
                for i = 3:length(columnName)

                    try
                        q.(columnName{i})(idx2) = NaN;
                    end
                    try
                        q.(columnName{i})(idx2) = NaT;
                    end
                    try
                        q.(columnName{i}){idx2} = {};
                    end
                end

                q.entryID{idx2} = x.entryID{s};
                q.participantID{idx2} = x.participantID{s};

                for i = 1:length(plist)

                    try
                        q.(plist{i})(idx2) = x.(plist{i})(s);
                    end
                    try
                        q.(plist{i})(idx2) = x.(plist{i})(s);
                    end
                    try
                        q.(plist{i}){idx2} = x.(plist{i}){s};
                    end

                end
            end
        end
    end
end
q = sortrows(q,"testTime_dsst","ascend");
writetable(q,fullfile('results','ect_01_octal.csv'));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = handleMultiEntry(x)
% rename subject based on their participant ID and timeStart
% x = alf_immediate;
sublist = unique(x.participantID);
y = table;
for s = 1:length(sublist)
    xx = x(find(strcmp(x.participantID,sublist{s})),:);
    if height(xx) == 1
        xx.entryID = {[sublist{s}]};
    else
        % lookfor testTime's column
        xx = sortrows(xx,find(contains(xx.Properties.VariableNames,"testTime")),"ascend");
        for ii = 1:height(xx)
            if ii == 1
                xx.entryID{ii} = [sublist{s}];
            else
                xx.entryID{ii} = [sublist{s} '_' num2str(ii)];
            end
        end
    end
    y = [y; xx];
end
y = movevars(y, "entryID", "Before", 1);
% alf_immediate = y;
end
