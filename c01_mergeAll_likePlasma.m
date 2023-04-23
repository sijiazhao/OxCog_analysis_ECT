clear; close all;
groupName = 'ProlificControl';

O = table;
D = readtable(fullfile('E:\GitHub\Masud_OxfordCognition\Project02_GiantProlific\demographics_oxcog_prolific','demographic_prolific.csv'));

% Merge OIS and OIS2
q = [readtable(fullfile('Results','result_OIS1.csv'));readtable(fullfile('Results','result_OIS2.csv'))];
writetable(q,fullfile('Results','result_OIS.csv'));


exptlist = {'REY','OIS','OMT','DSST','TMT','CORSI','TOL'};

% filename = 'merged_OIS';
for j = 1:length(exptlist)
    Q = readtable(fullfile('Results',['result_' exptlist{j} '.csv']));

    t = readtable(fullfile('ParameterNameList',['parameterList_task_' exptlist{j} '.xlsx']),'Sheet','all');
    plist = t.plist;
    typelist = t.typelist;

    if j == 1
        O.participantID = Q.participantID;
    end

    for i  = 1:length(plist)
        switch typelist{i}
            case 'num'
                O.(plist{i}) = nan(size(O.participantID));
            case 'time'
                %                 O.(plist{i}) = NaT(size(O.participantID));
                O.(plist{i}) = cell(size(O.participantID));
            case 'cell'
                O.(plist{i}) = cell(size(O.participantID));
        end
    end


    for s = 1:size(Q,1)
        idx = find(strcmp(O.participantID, Q.participantID{s}));
        if ~isempty(idx)
            for i = 1:length(plist)
                try
                    O.(plist{i})(idx) = Q.(plist{i})(s);
                catch
                    O.(plist{i}){idx} = Q.(plist{i}){s};
                end
            end

        else
            idx2 = size(O,1)+1;
            columnName = O.Properties.VariableNames;
            for i = 1:length(columnName)
%                 switch columnName{i}
%                     case {'testTime_OIS','testTime_OMT','testTime_CORSI','testTime_TMT','testTime_TOL','testTime_DSST','testTime_OG'}
%                         O.(columnName{i})(idx2) = NaT;
%                     otherwise
                        try
                            O.(columnName{i})(idx2) = NaN;
                        catch
                            O.(columnName{i}){idx2} = {};
                        end
%                 end
            end
            O.participantID{idx2} = Q.participantID{s};
            for i = 1:length(plist)
                try
                    O.(plist{i})(idx2) = Q.(plist{i})(s);
                catch
                    O.(plist{i}){idx2} = Q.(plist{i}){s};
                end
            end
        end
    end
    %     filename = [filename '_' exptlist{j}];
end
t = O;
for s = 1:height(t)
    idx = find(strcmp(D.participantID,t.participantID{s}));
    if ~isempty(idx)
        t.age(s) = D.age(idx);
        t.gender{s} = D.gender{idx};
    else
        t.age(s) = nan;
        t.gender{s} = 'Unknown';
    end
end
t = movevars(t,'gender','After','participantID');
t = movevars(t,'age','After','participantID');
O = t;

O.group = repmat({groupName},height(O),1);
O = movevars(O,'group','Before','participantID');
writetable(O,fullfile('all.csv'));
