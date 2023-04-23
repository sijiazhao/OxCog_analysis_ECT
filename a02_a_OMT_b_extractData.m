% After mixing 1 item and 3 items, I need to use different way to extract
% rows with n items
close all; clear;
addpath(genpath('Functions'));

chooseExperiment = 'OMTpavlovia';
repositoryName =  chooseExperiment;

path_in = fullfile('Data',repositoryName);
fileList = dir(fullfile(path_in,'*.csv'));
numSubjects = length(fileList);

R = table; % trial data
Q = table; % subject data

for s = 1:numSubjects
    filename = fullfile(fileList(s).folder,fileList(s).name);
    T = readtable(filename);

    experimentName = {'OMT_Pavlovia'};

    try
        participantID = T.('participantID')(1);
    catch
        participantID = T.('participantProlificPid')(1);
    end

    % Check IP
    userAgent = {T.participantBrowserUseragent{1}};

    % Check if the subject completed the experiment
    testTime = T.participantStartTime{1};
    testTime = strrep(testTime,'h','.');
    testTime = datetime(testTime,'InputFormat','yyyy-MM-dd_HH.mm.ss.SSS');
    testTime_OMT = testTime;
    taskCompletion = ~isempty(T.participantFinishTime{size(T,1)});

    if taskCompletion
        r = table; % trial info for this participant
        adjustRatio = T.adjustRatio(1);

        % Main 1item
        nItem = 1;
        r0 = extractData_OMTpavlovia(T,nItem);
        r = [r; r0];

        % Main 3item
        nItem = 3;
        r0 = extractData_OMTpavlovia(T,nItem);
        r = [r; r0];

        R = [R;r];

        % Summarise this participant
        ProportionCorrect_1item_4secs = nanmean(r.IdentificationResp(find(r.nItem == 1)));
        ProportionCorrect_3item_4secs = nanmean(r.IdentificationResp(find(r.nItem == 3)));

        AbsoluteError_1item_4secs = nanmean(r.AbsoluteError(find(r.nItem == 1 & r.IdentificationResp == 1)));
        AbsoluteError_3item_4secs = nanmean(r.AbsoluteError(find(r.nItem == 3 & r.IdentificationResp == 1)));

        %         IdentificationTime_1item_4secs = nanmean(fRemoveOutliers(r.IdentificationTime(find(r.nItem == 1))));
        %         IdentificationTime_3item_4secs = nanmean(fRemoveOutliers(r.IdentificationTime(find(r.nItem == 3))));
        %         LocalisationTime_1item_4secs = nanmean(fRemoveOutliers(r.LocalisationTime(find(r.nItem == 1))));
        %         LocalisationTime_3item_4secs = nanmean(fRemoveOutliers(r.LocalisationTime(find(r.nItem == 3))));

        %         IdentificationTime_1item_4secs = nanmean(r.IdentificationTime(find(r.nItem == 1 & r.IdentificationResp == 1)));
        %         IdentificationTime_3item_4secs = nanmean(r.IdentificationTime(find(r.nItem == 3 & r.IdentificationResp == 1)));
        %         LocalisationTime_1item_4secs = nanmean(r.LocalisationTime(find(r.nItem == 1 & r.IdentificationResp == 1)));
        %         LocalisationTime_3item_4secs = nanmean(r.LocalisationTime(find(r.nItem == 3 & r.IdentificationResp == 1)));

        IdentificationTime_1item_4secs = nanmean(fRemoveOutliers(r.IdentificationTime(find(r.nItem == 1 & r.IdentificationResp == 1))));
        IdentificationTime_3item_4secs = nanmean(fRemoveOutliers(r.IdentificationTime(find(r.nItem == 3 & r.IdentificationResp == 1))));
        LocalisationTime_1item_4secs = nanmean(fRemoveOutliers(r.LocalisationTime(find(r.nItem == 1 & r.IdentificationResp == 1))));
        LocalisationTime_3item_4secs = nanmean(fRemoveOutliers(r.LocalisationTime(find(r.nItem == 3 & r.IdentificationResp == 1))));

        % Compute misbinding rate using a permutation method (like OMT SERVER
        % SIMPLE)

        nTrial = size(r,1);
        samplenum = 5000;
        for i = 1:nTrial
            try
                misbindmatrix = [repmat(r.TargetDistance(i),samplenum,1) repmat(r.MisbindDistance(i),samplenum,1) datasample(rmmissing(r.MisbindDistance),samplenum)]; % [samplenum x 3] matrix of [targetDistance, nearest nonTarget distance, random trials' misbinding distances]
                [~,ii] = min(misbindmatrix,[],2); % which column is smallest
                r.misbinding(i)  = sum(ii==2)/samplenum; % mean rate of nonTargDistance being smallest
                r.guessrate(i)   = sum(ii==3)/samplenum; % mean rate of random trial nonTargDistance
                r.targetresp(i)  = sum(ii==1)/samplenum; % mean rate of target distance being smallest
                r.imprecision(i) = min(misbindmatrix(1,1:2)); % nearest neighbour distance
                clear misbindingmatrix ii
            catch % some subjects didn't localise at all.
                r.misbinding(i)  = nan; % mean rate of nonTargDistance being smallest
                r.guessrate(i)   = nan; % mean rate of random trial nonTargDistance
                r.targetresp(i)  = nan; % mean rate of target distance being smallest
                r.imprecision(i) = nan; % nearest neighbour distance
            end
        end
        Misbinding_Simple_1item_4secs = nanmean(r.misbinding(r.nItem==1));
        Misbinding_Simple_3item_4secs = nanmean(r.misbinding(r.nItem==3));
        Guessing_Simple_1item_4secs = nanmean(r.guessrate(r.nItem==1));
        Guessing_Simple_3item_4secs = nanmean(r.guessrate(r.nItem==3));
        Target_Simple_1item_4secs = nanmean(r.targetresp(find(r.nItem == 1)));
        Target_Simple_3item_4secs = nanmean(r.targetresp(find(r.nItem == 3)));
        Precision_Simu_1item_4secs = nan;
        Precision_Simu_3item_4secs = nan;
        Imprecision_Simple_1item_4secs = nanmean(r.imprecision(r.nItem==1));
        Imprecision_Simple_3item_4secs = nanmean(r.imprecision(r.nItem==3));

        userAgent_OMT = userAgent;

        screensize_OMT = r.screenWidth(1);

        q = table(participantID,userAgent_OMT, testTime_OMT, ...
            screensize_OMT,...
            ProportionCorrect_1item_4secs, ProportionCorrect_3item_4secs,AbsoluteError_1item_4secs,AbsoluteError_3item_4secs,...
            Target_Simple_1item_4secs,Target_Simple_3item_4secs,Misbinding_Simple_3item_4secs,Guessing_Simple_1item_4secs,Guessing_Simple_3item_4secs,...
            IdentificationTime_1item_4secs,IdentificationTime_3item_4secs,LocalisationTime_1item_4secs,LocalisationTime_3item_4secs,...
            Precision_Simu_1item_4secs,Precision_Simu_3item_4secs,Imprecision_Simple_1item_4secs,Imprecision_Simple_3item_4secs);
        Q = [Q; q];

    end
end
writetable(Q,fullfile('Results', 'result_OMTpavlovia.csv'));
mkdir(fullfile('Results','OMT'));
writetable(R,fullfile('Results','OMT','result_OMTpavlovia_trial.csv'));


% JOIN with server data
a = readtable(fullfile('Results','result_OMTserver.csv'));
v = Q.Properties.VariableNames;
v(strcmp(v,'screensize_OMT')) = [];
a = a(:,v);

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
        [~,i] = min(a.testTime_OMT);
        keepIdx = [keepIdx; idx(i)];

        idx(i) = [];
        removeIdx = [removeIdx; idx];
    end
end

if ~isempty(duplicates)
    duplicates

    Q(removeIdx,:) = [];
    disp('Duplicates have been removed. Only first complete entry kept.');
else
    disp('No duplicates!')
end

Q.OMT_ProportionCorrect = mean([Q.ProportionCorrect_1item_4secs,Q.ProportionCorrect_3item_4secs],2);
Q.OMT_AbsoluteError = mean([Q.AbsoluteError_1item_4secs,Q.AbsoluteError_3item_4secs],2);
Q.OMT_TargetDetection = mean([Q.Target_Simple_1item_4secs,Q.Target_Simple_3item_4secs],2);
Q.OMT_Misbinding = Q.Misbinding_Simple_3item_4secs;
Q.OMT_Guessing = mean([Q.Guessing_Simple_1item_4secs,Q.Guessing_Simple_3item_4secs],2);
Q.OMT_Precision = mean([Q.Precision_Simu_1item_4secs,Q.Precision_Simu_3item_4secs],2);
Q.OMT_Imprecision = mean([Q.Imprecision_Simple_1item_4secs,Q.Imprecision_Simple_3item_4secs],2);
Q.OMT_IdentificationTime = mean([Q.IdentificationTime_1item_4secs,Q.IdentificationTime_3item_4secs],2);
Q.OMT_LocalisationTime = mean([Q.LocalisationTime_1item_4secs,Q.LocalisationTime_3item_4secs],2);

writetable(Q,fullfile('Results', 'result_OMT.csv'));






%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = fRemoveOutliers(a)
% a(isnan(a)) = [];
a(a>nanmean(a)+2*nanstd(a) | a<nanmean(a) - 2*nanstd(a)) = nan;
end


function r0 = extractData_OMTpavlovia(T,nItem)
r0 = table;
screenUnit = 'px'; % raw, cm, px

% Get screensize from card calibration (if applicable
pixels = T.pixels;
pixels = pixels(~cellfun(@isempty,pixels)); pixels = pixels{1};

a = jsondecode(pixels);
pxWidth = a.x_pxWidth;
pxHeight = a.x_pxHeight;

cmWidth = a.x_cmWidth;
cmHeight = a.x_cmHeight;

screenWidth_cm = cmWidth;
screenHeight_cm = cmHeight;

switch screenUnit
    case 'px'
        screenWidth = pxWidth;
        screenHeight = pxHeight;
    case 'cm'
        screenWidth = cmWidth;
        screenHeight = cmHeight;
end

% Localise the rows for trials
tmp = T.experiment;
tmp = tmp(~cellfun(@isempty,tmp));
tmp = tmp{1};
experimentVersion = tmp;
switch experimentVersion
    case 'OMT' % experiment = 'OMT'
        % Early version used block design: 1 item 1 lock, 3 item 1 block.
        % so each block has a loop name containing item number
        pname = ['loop_main_' num2str(nItem) 'item_thisTrialN'];
        t1 = T(~isnan(T.(pname)),:);
        t1.trialIndex = t1.(pname)+1;
        t1.nItem = ones(size(t1,1),1) * nItem;
    otherwise
        t1 = T(~isnan(T.('loop_main_3item_thisTrialN')),:);
        t1.trialIndex = t1.('loop_main_3item_thisTrialN')+1;
        t1.nItem = nan(size(t1,1),1);
        for ii = 1:size(t1,1)
            t1.nItem(ii) = ~isnan(t1.item1_shapeIdx(ii)) + ~isnan(t1.item2_shapeIdx(ii)) + ~isnan(t1.item3_shapeIdx(ii));
        end

        t1 = t1(t1.nItem == nItem,:);
end

try
    participantID = T.participantName(1);
catch
    participantID = T.participantProlificPid(1);
end

trialIdxList = unique(t1.trialIndex)';
for i = trialIdxList
    a = t1(t1.trialIndex == i,:);

    % find the last entry where the column object contains the string
    a(strcmp(a.object,''),:) = [];
    ai = a(size(a,1),:); % Last row contains most information of this trial

    c = table;
    c.participantID = participantID;
    c.trialIdx = i;
    c.nItem = nItem;

    pname = {'item1_locationX','item1_locationY','item2_locationX','item2_locationY','item3_locationX','item3_locationY',...
        'item1','item2','item3','target','target_itemIdx'};
    for ii = 1:length(pname)
        switch pname{ii}
            case {'item2','item3'}
                b = a.(pname{ii})(1);
                try
                    if isnan(b)
                        c.(pname{ii}) = {''};
                    end
                catch
                    c.(pname{ii}) = a.(pname{ii})(1);
                end
            otherwise
                c.(pname{ii}) = a.(pname{ii})(1);
        end
    end


    c.IdentificationTime = a.time(1);
    c.LocalisationTime = ai.time;
    c.objTarget = ai.target{1};
    c.objResponse = ai.object{1};
    % Chose the right object?
    if strcmp(c.objResponse, c.objTarget)
        c.IdentificationResp = 1;
    else
        c.IdentificationResp = 0;
    end

    % Get coordinates for object and response
    try
        c.finalLocX = ai.objectX;
        c.finalLocY = ai.objectY;
    catch
        a = t1(t1.trialIndex == i,:);
        a(isnan(a.x),:) = [];
        ai.x = a.x(size(a,1));
        ai.y = a.y(size(a,1));
        c.finalLocX = ai.x;
        c.finalLocY = ai.y;
    end
    c.locTargetX = c.(['item' num2str(c.target_itemIdx) '_locationX']);
    c.locTargetY = c.(['item' num2str(c.target_itemIdx) '_locationY']);

    switch screenUnit
        case {'px','cm'}
            c.finalLocX = c.finalLocX*screenWidth;
            c.finalLocY = c.finalLocY*screenHeight;
            c.locTargetX = c.locTargetX*screenWidth;
            c.locTargetY = c.locTargetY*screenHeight;
    end
    p1 = [c.locTargetX,c.locTargetY]; p2 = [c.finalLocX,c.finalLocY];
    er = norm(p1-p2);
    c.AbsoluteError = er;

    if nItem == 3

        % Option B: Misbinding object and location
        p1 = [c.finalLocX,c.finalLocY]; % response location
        p2 = [...
            c.item1_locationX,c.item1_locationY;...
            c.item2_locationX,c.item2_locationY;...
            c.item3_locationX,c.item3_locationY;...
            ];

        switch screenUnit
            case {'px','cm'}
                p2(:,1) = p2(:,1)*screenWidth;
                p2(:,2) = p2(:,2)*screenHeight;
        end
        dist = nan(nItem,1);
        for jj = 1:nItem
            dist(jj) = pdist([p1; p2(jj,:)]);
        end
        ik = c.target_itemIdx; % if successfully binding, this index should be the minimum distance


        c.TargetDistance = dist(ik);

        dist(ik) = [];
        [minDist, minIdx] = min(dist);


        c.MisbindDistance = minDist; % smallest distance to nonTarget
    else
        c.TargetDistance = c.AbsoluteError;
        c.MisbindDistance = nan;
    end

    c.screenWidth = screenWidth_cm;
    c.screenHeight = screenHeight_cm;

    r0 = [r0; c];
end
end
