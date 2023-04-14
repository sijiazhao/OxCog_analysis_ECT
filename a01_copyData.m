clear; close all;
T = readtable(fullfile('copyData.xlsx'));

currentFolder = pwd;
[~,currentFolderName,~] = fileparts(pwd);
path_root = strrep(pwd,currentFolderName,''); % get the path of the parent folder

% for i = 2:size(T,1)
for i = 1:size(T,1)
    path_out = fullfile('Data',T.path_out{i});
    try
        rmdir(path_out);
    end
    mkdir(path_out);

    nameProject = T.project{i};
    nameTask = T.task{i};

    path_in = fullfile(path_root,T.path_in{i},'data'); % where is the data?
    fileList = dir(fullfile(path_in,'*.csv'));
    n = length(fileList);
    for s = 1:n
        filename = fullfile(fileList(s).folder,fileList(s).name);

        switch nameTask
            case {'TMT','TOL','DSST'}
                tbl = readtable(filename, 'Delimiter', ',');
            otherwise
                tbl = readtable(filename);
        end

        if height(tbl) > 2 % nothing in the file

            switch nameTask
                case 'TOL'
                    taskCompletion = sum(strcmp(tbl.event,'task finish'));
                otherwise
                    taskCompletion = ~isempty(tbl.participantFinishTime{size(tbl,1)});
            end

            if taskCompletion % copy if complete
                % --------------------------------------------------------------------%
                % Correct Participant ID
                participantID = tbl.('participantName'){1};


                participantID = strrep(participantID,'.','_');
                participantID = strrep(participantID,'-','_');

                participantID = strrep(participantID,'plasma','');
                participantID = strrep(participantID,'Plasma','');
                participantID = strrep(participantID,'ST_','');
                participantID = strrep(participantID,'st_','');
                participantID = strrep(participantID,'ST','');
                participantID = strrep(participantID,'st','');

                participantID = strrep(participantID,'\','');
                participantID = strrep(participantID,'/','');

                participantID = strrep(participantID,'p','P');
                participantID = strrep(participantID,'c','C');
                participantID = strrep(participantID,'v','');
                participantID = strrep(participantID,'V','');
                participantID = strrep(participantID,' ','');
                participantID = strrep(participantID,'__','_');

                if strcmp(participantID(length(participantID)),'_')
                    participantID(length(participantID)) = [];
                end

                tbl.('participantID'){1} = participantID;
                % --------------------------------------------------------------------%

                try
                    testDate = tbl.participantStartTime{1};
                catch % TOL
                    testDate = tbl.time{1};
                end

                filename2 = [participantID,'_', testDate, '.csv'];
                writetable(tbl,fullfile(path_out,filename2));

                %         copyfile(filename, new_path_in);
                disp(filename2);
            end
        end
    end
end
