% Version 20220817 - Remove OMT misbinding, keep visual misbinding =
% identify a familiar object to a wrong scene.

% Version 20220412 - Different ways to compute misbinding rate
% 1. Old way: choose wrong object which is an object that you have seen
% before
% 2. OMT way: misbinding the location and object: wrong object but right
% position.

close all; clear;
addpath(genpath('Functions'));
mkdir('Results');

sublist_completeStatus = table;

ccccount = 0;

chooseExperiment = 'OIS';
nTrial = 20;
nBlock = 4;
repositoryName =  chooseExperiment;
distanceUnit = 'cm'; % 'pixel','percent'


path_in = fullfile('Data',repositoryName);
conditionTable = readtable(fullfile('Data','condition_OIS.xlsx'));
for i = 1:size(conditionTable,1)
    eval(['loc = ' conditionTable.loc{i} ';']);
    conditionTable.x(i) = loc(1);
    conditionTable.y(i) = loc(2);
end
mapme = readtable(fullfile('Data','mapObject_OIS.xlsx'));
if size(mapme,1) ~= nTrial; error('Wrong mapObject.xlsx!'); return; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Extract trial data
fileList = dir(fullfile(path_in,'*.csv'));
numSubjects = length(fileList);

R = table; % single trial
Q = table; % one row for each subject
B = table; % one row for one block for each subject

n = 0;
for s = 1:numSubjects
    filename = fullfile(fileList(s).folder,fileList(s).name);
    T = readtable(filename);

    if size(T,1) == 1 && size(T,2) > 200
        T = readtable(filename, 'Delimiter', ',');
    end


    try
        opts = detectImportOptions(filename);
        T2 = readtable(filename,opts);
        T.realtime = T2.realtime;
        T = movevars(T,'realtime','Before','Var1');
        n = n+1;
    end


    experimentName = {chooseExperiment};

    try
        participantID = T.('participantID')(1);
    catch
        participantID = T.('participantProlificPid')(1);
    end

    % Check IP
    userAgent = {T.participantBrowserUseragent{1}};

    % Time start the task
    testTime = T.participantStartTime{1};
    testTime = strrep(testTime,'h','.');
    testTime = datetime(testTime,'InputFormat','yyyy-MM-dd_HH.mm.ss.SSS');

    % If participant completed the task
    taskCompletion = ~isempty(T.participantFinishTime{size(T,1)});

    if taskCompletion

        % Check screensize and result from card calibration
        if isempty(T.pixels{1})
            a = jsondecode(T.pixels{2});
        else
            a = jsondecode(T.pixels{1});
        end
        try
            pxWidth = a.x_pxWidth;
            pxHeight = a.x_pxHeight;
            cmWidth = a.x_cmWidth;
            cmHeight = a.x_cmHeight;
        catch
            pxWidth = a.display_i.width;
            pxHeight = a.display_i.height;
            cmWidth = a.display_cm.width;
            cmHeight = a.display_cm.height;
        end

        % Compute the actual locations based on the condition table
        % (Yes, I also saved the actual locations in the results,
        % but here because I need to access to all in this block so
        % I need to do in this way)
        conditionTable2 = conditionTable;
        conditionTable2.x2 = nan(size(conditionTable,1),1);
        conditionTable2.y2 = nan(size(conditionTable,1),1);
        adjustRatio = T.adjustRatio(1);
        offset = 0; % check psychopy code (in LoadJSON, setOffset component)

        conditionTable2.x2 = conditionTable2.x*adjustRatio;
        conditionTable2.y2 = conditionTable2.y*adjustRatio + offset;

        %% Trial analysis - the phase of ImmediateRecall
        phaseName = 'ImmediateRecall';
        for block = 1:nBlock
            for i = 0:(nTrial/nBlock-1)
                extra_OIS_analyseSingleTrial;
            end
        end

        %%  Trial analysis - DelayedRecall
        phaseName = 'DelayedRecall';
        for block = 1:nBlock
            for i = 0:(nTrial/nBlock-1)
                extra_OIS_analyseSingleTrial;
            end
        end

        %% Block analysis
        for block = 1:nBlock
            idx = find(strcmp(R.participantID, participantID{1}) & strcmp(R.phase,'ImmediateRecall') & R.block == block & R.validTrial == 1);
            ImmediateObjectAccuracy = nanmean(R.DelayedObjectAccuracy(idx));
            ImmediateSemanticAccuracy = nanmean(R.DelayedSemanticAccuracy(idx));
            ImmediateLocationError = nanmean(R.LocationError(idx));

            idx = find(strcmp(R.participantID, participantID{1}) & strcmp(R.phase,'DelayedRecall') & R.block == block & R.validTrial == 1);
            DelayedObjectAccuracy = nanmean(R.DelayedObjectAccuracy(idx));
            DelayedSemanticAccuracy = nanmean(R.DelayedSemanticAccuracy(idx));
            DelayedLocationError = nanmean(R.LocationError(idx));

            Mismatching_FamiliarObject = nanmean(R.Mismatching_FamiliarObject(idx));
            Guessing_NovelObject = nanmean(R.Guessing_NovelObject(idx));
            Misbinding_ObjectLocation = nanmean(R.Misbinding_ObjectLocation(idx));
            Guessing_WrongObjectRandomLocation = nanmean(R.Guessing_WrongObjectRandomLocation(idx));

            TargetDetection = nanmean(R.TargetDetection(idx));

            u = symunit;
            cm2inch = double(unitConversionFactor(u.cm,u.inch));
            screensize = cm2inch * sqrt(T.display_cm_height(1)^2 + T.display_cm_width(1)^2);

            B = [B; table(block,experimentName,participantID,...
                ImmediateObjectAccuracy, DelayedObjectAccuracy,ImmediateLocationError,DelayedLocationError, ...
                Mismatching_FamiliarObject, Guessing_NovelObject,TargetDetection,...
                ImmediateSemanticAccuracy,DelayedSemanticAccuracy, screensize)];
        end

        %         idx = find(strcmp(R.participantID, participantID{1}) & strcmp(R.phase,'ImmediateRecall'));
        idx = find(strcmp(R.participantID, participantID{1}) & strcmp(R.phase,'ImmediateRecall') & R.validTrial == 1);
        ImmediateObjectAccuracy = nanmean(R.DelayedObjectAccuracy(idx))*100;
        ImmediateSemanticAccuracy = nanmean(R.DelayedSemanticAccuracy(idx))*100;
        ImmediateLocationError = nanmean(R.LocationError(idx));


        %         idx = find(strcmp(R.participantID, participantID{1}) & strcmp(R.phase,'DelayedRecall'));
        idx = find(strcmp(R.participantID, participantID{1}) & strcmp(R.phase,'DelayedRecall') & R.validTrial == 1);
        DelayedObjectAccuracy = nanmean(R.DelayedObjectAccuracy(idx))*100;
        DelayedSemanticAccuracy = nanmean(R.DelayedSemanticAccuracy(idx))*100;
        DelayedLocationError = nanmean(R.LocationError(idx));

        Mismatching_FamiliarObject = nanmean(R.Mismatching_FamiliarObject(idx))*100;
        Guessing_NovelObject = nanmean(R.Guessing_NovelObject(idx))*100;
        Misbinding_ObjectLocation = nanmean(R.Misbinding_ObjectLocation(idx))*100;
        Guessing_WrongObjectRandomLocation = nanmean(R.Guessing_WrongObjectRandomLocation(idx))*100;
        TargetDetection = nanmean(R.TargetDetection(idx))*100;

        u = symunit;
        cm2inch = double(unitConversionFactor(u.cm,u.inch));
        screensize_OIS = cm2inch * sqrt(T.display_cm_height(1)^2 + T.display_cm_width(1)^2);

        testTime_OIS = testTime;
        userAgent_OIS = userAgent;
    else
        ImmediateObjectAccuracy = nan;
        DelayedObjectAccuracy = nan;
        Mismatching_FamiliarObject = nan;
        Guessing_NovelObject = nan;
        Misbinding_ObjectLocation = nan;
        Guessing_WrongObjectRandomLocation = nan;
        TargetDetection = nan;
        ImmediateSemanticAccuracy = nan;
        DelayedSemanticAccuracy = nan;
        ImmediateLocationError = nan;
        DelayedLocationError = nan;
        screensize_OIS = nan;

        testTime_OIS = testTime;
        userAgent_OIS = userAgent;
    end

    Q = [Q; table(participantID, testTime_OIS, userAgent_OIS, screensize_OIS,...
        ImmediateObjectAccuracy,DelayedObjectAccuracy,ImmediateLocationError,DelayedLocationError,...
        Mismatching_FamiliarObject, Guessing_NovelObject, TargetDetection, ...
        Misbinding_ObjectLocation, Guessing_WrongObjectRandomLocation,...
        ImmediateSemanticAccuracy,DelayedSemanticAccuracy)];
end

Q.change_in_ObjectAccuracy = -(Q.ImmediateObjectAccuracy - Q.DelayedObjectAccuracy);
Q.change_in_SemanticAccuracy = -(Q.ImmediateSemanticAccuracy - Q.DelayedSemanticAccuracy);
Q.change_in_LocationError = -(Q.ImmediateLocationError - Q.DelayedLocationError);


%% Remove the participant who didn't learn (accuracy for correct object < 20%)
thresholdLearning = 0.5;
badsub = find(Q.ImmediateObjectAccuracy < thresholdLearning);
Qbad = Q(badsub,:);
% Q(badsub,:) = [];


filename = fullfile('Results','result_OIS.csv');
writetable(Q,filename);

mkdir(fullfile('Results','OIS'));
filename = fullfile('Results','OIS','individual.csv');
writetable(R,filename);
if nBlock >1
    filename = fullfile('Results','OIS','block.csv');
    writetable(B,filename);
end


% More checks like duplicates, task completion etc.
% Check if there are any duplicates
sublist = Q.participantID;
duplicates = table;
uniqsublist  = unique(sublist);
keepIdx = [];
removeIdx = [];
for s = 1:length(uniqsublist)
    idx = find(strcmp(sublist,uniqsublist{s}));
    if length(idx) == 1
        keepIdx = [keepIdx; idx];
    elseif length(idx) > 1
        subject = {uniqsublist{s}};
        duplicates = [duplicates; table(subject)];
        a = Q(idx,:);

        % Keep the first complete entry
        [~,i] = min(a.testTime_OIS);
        keepIdx = [keepIdx; idx(i)];

        idx(i) = [];
        removeIdx = [removeIdx; idx];
    end
end

if ~isempty(duplicates)
    duplicates
end
Q(removeIdx,:) = [];
disp('Duplicates have been removed. Only first complete entry kept.');

filename = fullfile('Results','result_OIS.csv');


% Rename parameters
% Q = readtable(fullfile('Results','result_OIS.csv'));
oldnames = {'ImmediateObjectAccuracy','DelayedObjectAccuracy','change_in_ObjectAccuracy',...
    'ImmediateLocationError','DelayedLocationError','change_in_LocationError',...
    'ImmediateSemanticAccuracy','DelayedSemanticAccuracy','change_in_SemanticAccuracy',...
    'Mismatching_FamiliarObject', 'Guessing_NovelObject', 'Misbinding_ObjectLocation', 'Guessing_WrongObjectRandomLocation','TargetDetection'};
newnames = {};
for i = 1:length(oldnames)
    newnames{i} = ['OIS_' oldnames{i}];
end
Q = renamevars(Q, oldnames, newnames);

% X = [];
% plist = {'OIS_ImmediateObjectAccuracy','OIS_DelayedObjectAccuracy','OIS_TargetDetection'};
% for i = 1:length(plist)
%     X = [X, Q.(plist{i})];
% end
% 
% [coeff,score,latent] = pca(X);
% Q.OIS_PCA = score(:,1);


% a = Q.OIS_PCA;
% b = Q.OIS_DelayedLocationError;
% ri = isnan(a) | isnan(b);
% a(ri) = []; b(ri) = [];
% [r,p] = corr(a,b)

writetable(Q,fullfile('Results','result_OIS.csv'));



% mkdir('Status');
% filename = fullfile('Status','completeStatus_OIS.csv');
% writetable(sublist_completeStatus,filename);


figure(1);
subplotx = 2;
subploty = 2;
i = 1;
subplot(subplotx,subploty,i);
histogram(Q.OIS_ImmediateObjectAccuracy);
title('Immediate Object Accuracy');

i = i+1;
subplot(subplotx,subploty,i);
histogram(Q.OIS_DelayedObjectAccuracy);
title('Delayed Object Accuracy');

i = i+1;
subplot(subplotx,subploty,i);
histogram(Q.OIS_ImmediateLocationError);
title('Immediate Location Error');

i = i+1;
subplot(subplotx,subploty,i);
histogram(Q.OIS_DelayedLocationError);
title('Delayed Location Error');
