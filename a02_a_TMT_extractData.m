close all; clear;
addpath(genpath('Functions'));
chooseExperiment = 'TMT';
path_in = fullfile('Data',chooseExperiment);
mkdir(fullfile('Results',chooseExperiment));

fileList = dir(fullfile(path_in,'*.csv'));
numSubjects = length(fileList);

R = table;
Q = table;

nTrial = 3;     % Part A: 1,2,3,4... % part = 1, iTrial = 1,2,3

for s = 1:numSubjects
    filename = fullfile(fileList(s).folder,fileList(s).name);
    T = readtable(filename, 'Delimiter', ',');

    subject = T.('participantID'){1};
    participantID = {subject};

    testTime = T.participantStartTime{1};
    testTime = strrep(testTime,'h','.');
    testTime = datetime(testTime,'InputFormat','yyyy-MM-dd_HH.mm.ss.SSS');
    testTime_TMT = testTime;

    taskCompletion = ~isempty(T.participantFinishTime{size(T,1)});
    if taskCompletion

        % Control: Connect two diagnal circles; 4 trials
        for i = 1:4
            idx = find(T.part == 0 & T.iTrial == i);
            trail = {'Control'};
            trialIdx = i;
            timeTaken = T.time(idx);
            R = [R; table(participantID, trail, trialIdx, timeTaken)];
        end

        % Part A: 1,2,3,4... % part = 1, iTrial = 1,2,3
        for i = 1:nTrial
            idx = find(T.part == 1 & T.iTrial == i);
            trail = {'Trail A'};
            trialIdx = i;
            timeTaken = T.time(idx);
            R = [R; table(participantID, trail, trialIdx, timeTaken)];
        end

        % Part B: 1,A,2,B,3,C, .. % part = 2, iTrial = 1,2,3
        for i = 1:nTrial
            idx = find(T.part == 2 & T.iTrial == i);
            trail = {'Trail B'};
            trialIdx = i;
            timeTaken = T.time(idx);
            R = [R; table(participantID, trail, trialIdx, timeTaken)];
        end

        c = R.timeTaken(strcmp(R.participantID, subject) & strcmp(R.trail,'Control'));
        a = R.timeTaken(strcmp(R.participantID, subject) & strcmp(R.trail,'Trail A'));
        b = R.timeTaken(strcmp(R.participantID, subject) & strcmp(R.trail,'Trail B'));

        c(c>mean(c)+2*std(c)) = [];
        a(a>mean(a)+2*std(a)) = [];
        b(b>mean(b)+2*std(b)) = [];

        TMT_A = mean(a);
        TMT_B = mean(b);
        TMT_control = mean(c);
        TMT_average =mean([TMT_A, TMT_B]);
        TMT_BdA = TMT_B/TMT_A;
        Q = [Q; table(participantID, testTime_TMT, TMT_A, TMT_B, TMT_control, TMT_average,TMT_BdA)];
    end
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
        [~,i] = min(a.testTime_TMT);
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


% Q.TMT_A_divide_control = Q.TMT_A./Q.TMT_control;
% Q.TMT_B_divide_control = Q.TMT_B./Q.TMT_control;

filename = fullfile('Results',['result_' chooseExperiment '.csv']);
writetable(Q,filename);
filename = fullfile('Results',chooseExperiment,['result_' chooseExperiment '_individual.csv']);
writetable(R,filename);
