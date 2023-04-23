% analyse DSST
% Created on 20220710

close all; clear;
addpath(genpath('Functions'));

chooseExperiment = 'DSST';

repositoryName = chooseExperiment;
path_in = fullfile('Data',repositoryName);

fileList = dir(fullfile(path_in,'*.csv'));
numSubjects = length(fileList);


R = table;
Q = table;

RT = {}; % PLOT LEARNING RATE: RESPONSE TIME OVER 120 S
cRT = {}; % PLOT LEARNING RATE: CORRECT RESPONSE TIME OVER 120 S
for s = 1:numSubjects
    filename = fullfile(fileList(s).folder,fileList(s).name);
    T = readtable(filename, 'Delimiter', ',');
    experimentName = {chooseExperiment};

    testTime = T.participantStartTime{1};
    testTime = strrep(testTime,'h','.');
    testTime_DSST = datetime(testTime,'InputFormat','yyyy-MM-dd_HH.mm.ss.SSS');
    taskCompletion = ~isempty(T.participantFinishTime{size(T,1)});


    if taskCompletion

        subject = T.('participantID'){1};


        idx = find(~isnan(T.responseCount)); idx = idx(2);
        T2 = T(idx,:);


        symbols = T2.symbols{1};
        symbols = strrep(symbols,'(','[');
        symbols = eval(strrep(symbols,')',']'));
        symbols = symbols';
        targets = T2.targets{1}; targets = strrep(targets,'(','{'); targets = eval(strrep(targets,')','}')); targets = targets';
        responses = eval(T2.responses{1}); responses = responses';
        responseTimes = eval(T2.responseTimes{1}); responseTimes = responseTimes';
        rt = [0; responseTimes]; rt = diff(rt);

        targets = targets(1:length(responses));

        index = nan(size(targets));
        for i = 1:length(targets)
            index(i) = find(strcmp(symbols,targets{i}));
        end

        participantID = repmat({subject},size(targets));
        time = repmat(testTime,size(targets));
        r = table(participantID,time,targets,index,responses,rt);
        R = [R; r];




        DSST_nAllResponse = T2.responseCount;
        DSST_nCorrectResponse = sum(index == responses);
        DSST_pHit = mean(index == responses);

        rt = rt(index == responses);
        rt(rt > mean(rt)+2*std(rt) | rt < mean(rt)-2*std(rt)) = [];
        DSST_rt = mean(rt);

        participantID = {subject};

        Q = [Q; table(participantID,testTime_DSST, DSST_nCorrectResponse,DSST_nAllResponse,DSST_pHit,DSST_rt)];

        RT = [RT; {responseTimes}];
        cRT = [cRT; {responseTimes(index == responses)}];
    end
end


% Check duplicates and only keep the first entry
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
        [~,i] = min(a.testTime_DSST);
        keepIdx = [keepIdx; idx(i)];

        idx(i) = [];
        removeIdx = [removeIdx; idx];
    end
end

if ~isempty(duplicates)
    duplicates
end

Q(removeIdx,:) = [];
RT(removeIdx,:) = [];
cRT(removeIdx,:) = [];
disp('Duplicates have been removed. Only first complete entry kept.');


%% Plot summary of DSST
figure(1);
fx = 2;
fy = 3;

i = 1;
subplot(fx,fy,i);
histogram(Q.DSST_nCorrectResponse);
title('Number of correct responses');
axis square;

i = i+1;
subplot(fx,fy,i);
histogram(Q.DSST_pHit*100);
title('Hit rate (%)');
axis square;

i = i+1;
subplot(fx,fy,i);
histogram(Q.DSST_rt);
title('Mean RT (s)');
axis square;

subplot(2,1,2); hold on;
for i = 1:length(RT)
    x = RT{i};
    y = i*ones(size(x));
    s = scatter(x,y,'.');
    s.MarkerFaceColor = 'k';
    s.MarkerEdgeColor = 'k';

    x1 = cRT{i};
    y1 = i*ones(size(x1));
    s = scatter(x1,y1,'.');
    s.MarkerFaceColor = 'g';
    s.MarkerEdgeColor = 'g';
end
xlim([0,121]);
ylim([0,i+1]);
xlabel('Time (s)');
ylabel('Participant index');

writetable(Q,fullfile('Results',['result_' chooseExperiment '.csv']));

