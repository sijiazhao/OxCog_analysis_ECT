clear; close all;

% O = readtable(fullfile('Results','result_OIS.csv'));O = O(:,1:11);

O = table;
D = readtable(fullfile('demographic.csv'));

exptlist = {'REY','OIS','OMT','DSST','TMT','CORSI','TOL'};

% filename = 'merged_OIS';
for j = 1:length(exptlist)
    Q = readtable(fullfile('Results',['result_' exptlist{j} '.csv']));

    t = readtable(fullfile('ParameterNameList',['parameterList_task_' exptlist{j} '.xlsx']),'Sheet','all');
    plist = t.plist;
    typelist = t.typelist;


    % Correct participant ID
    t = readtable('E:\Dropbox (Neurological Conditions)\Plasma Biomarkers\sijiazhao_correctParticipantID_plasma.xlsx');
    for s = 1:height(Q)
        idx = find(strcmp(t.old, Q.participantID{s}));
        if ~isempty(idx)
            disp([Q.participantID{s} ' --> ' t.new{idx}]);
            Q.participantID{s} = t.new{idx};
        end
    end

    t = readtable('E:\Dropbox (Neurological Conditions)\Plasma Biomarkers\sijiazhao_rejectParticipantID_plasma.xlsx');
    rejectList = [];
    for s = 1:height(Q)
        idx = find(strcmp(t.reject, Q.participantID{s}));
        if ~isempty(idx)
            rejectList = [rejectList; s];
        end
    end
    Q(rejectList,:) = [];

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
        t.group{s} = D.group{idx};
    else
        t.age(s) = nan;
        t.gender{s} = 'unknown';
        %         switch t.participantID{s}(1)
        %             case 'P'
        %                 t.group{s} = 'AD';
        %             case 'C'
        %                 t.group{s} = 'Control';
        %             otherwise
        t.group{s} = 'unknown';
        %         end
    end
end
t = movevars(t,'gender','After','participantID');
t = movevars(t,'age','After','participantID');
O = t;



for s = 1:height(O)
    subject = O.participantID{s};
    newstr = split(subject,'_');
    O.cc_code{s} = newstr{1};
    try
        O.visit(s) = newstr{2};
    catch
        O.visit(s) = nan;
    end
end
O = movevars(O,'group','Before','participantID');
O = movevars(O,'visit','Before','participantID');
O = movevars(O,'cc_code','After','participantID');
O = sortrows(O,'visit','ascend');
O = sortrows(O,"group","ascend");
writetable(O,fullfile('OnlineTaskCompletionStatus_PlasmaProject.xlsx'));
