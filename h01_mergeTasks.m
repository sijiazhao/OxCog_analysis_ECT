close all; clear;
try
    addpath(genpath('E:\GitHub\Functions'));
end

path_root = 'E:\GitHub\Masud_OxfordCognition\Project13_ECT';

alf_immediate = readtable(fullfile(path_root, 'ect_alf', 'results', 'a.csv')); alf_immediate = handleMultiEntry(alf_immediate);
alf_delayed = readtable(fullfile(path_root, 'ect_alf_delayed', 'results', 'a.csv')); alf_delayed = handleMultiEntry(alf_delayed);

ois = readtable(fullfile(path_root, 'ect_ois', 'results', 'a.csv')); ois = handleMultiEntry(ois);
omt = readtable(fullfile(path_root, 'ect_omt', 'results', 'a.csv')); omt = handleMultiEntry(omt);

tmt = readtable(fullfile(path_root, 'ect_tmt', 'results', 'a.csv')); tmt = handleMultiEntry(tmt);
dsst = readtable(fullfile(path_root, 'ect_dsst', 'results', 'a.csv')); dsst = handleMultiEntry(dsst);
srt = readtable(fullfile(path_root, 'ect_srt', 'results', 'a.csv')); srt = handleMultiEntry(srt);
crt = readtable(fullfile(path_root, 'ect_crt', 'results', 'a.csv')); crt = handleMultiEntry(crt);
gonogo = readtable(fullfile(path_root, 'ect_gonogo', 'results', 'a.csv')); gonogo = handleMultiEntry(gonogo);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% merge all tables
taskList = {'dsst','alf_immediate','alf_delayed','ois','omt','tmt','srt','crt','gonogo'};
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
            q.participantID_unique = x.participantID_unique;
                    q.participantID = x.participantID;
        end

        for i  = 1:length(plist)
            switch typelist{i}
                case 'num'
                    q.(plist{i}) = nan(size(q.participantID_unique));
                case 'time'
                    q.(plist{i}) = NaT(size(q.participantID_unique));
                case 'cell'
                    q.(plist{i}) = cell(size(q.participantID_unique));
            end
        end

        for s = 1:size(x,1)
            idx = find(strcmp(q.participantID_unique, x.participantID_unique{s}));
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

                q.participantID_unique{idx2} = x.participantID_unique{s};
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
writetable(q,'OxCog_ECT.csv');

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
        xx.participantID_unique = {[sublist{s} '_1']};
    else
        % lookfor testTime's column
        xx = sortrows(xx,find(contains(xx.Properties.VariableNames,"testTime")),"ascend");
        for ii = 1:height(xx)
            xx.participantID_unique{ii} = [sublist{s} '_' num2str(ii)];
        end
    end
    y = [y; xx];
end
y = movevars(y, "participantID_unique", "Before", 1);
% alf_immediate = y;
end
