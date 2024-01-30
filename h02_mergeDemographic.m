clear;
d = readtable(fullfile('ParticipantInformation','DNS_Battery_Batchfile.xlsx'));
q = readtable(fullfile('results','ect_01_octal.csv'));


% match q.entryID with d.entryID
q.weekID = nan(height(q),1);
q.participant = cell(height(q),1);
for s = 1:height(q)
    idx = find(strcmp(d.entryID,q.entryID{s}));
    if ~isempty(idx)
        q.participant{s} = d.participantID{idx};
        q.weekID(s) = d.weekID(idx);
    end
end

d.dataExist_any = zeros(height(d),1);
d.dataExist_full_1 = zeros(height(d),1);
d.dataExist_full_2 = zeros(height(d),1);
for s = 1:height(d)
    idx = find(strcmp(q.entryID,d.entryID{s}));
    if ~isempty(idx)
        d.dataExist_any(s) = 1;
    end
    if ~isempty(idx) && ~isnan(q.ALF_nDelayedRecall(idx))
        d.dataExist_full_1(s) = 1;
    end
    if ~isempty(idx) && ~isnan(q.srt_RT(idx))
        d.dataExist_full_2(s) = 1;
    end
end
q = movevars(q, "weekID", "Before", "entryID");
q = movevars(q, "participant", "Before", "weekID");
q = movevars(q, "entryID", "Before", "participant");

writetable(q,fullfile('results','ect_02a_octal_participantID_matched.csv'));
writetable(d,fullfile('results','ect_02b_demo_participantID_matched.csv'));

