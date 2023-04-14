function OMTFRACTALS1AND3
[yourdirectory,~,~] = fileparts(mfilename('fullpath')); %find out where this function is
root  = yourdirectory; %make it your root
addpath(genpath(root)) %add all subfolders to matlab path
cd(root) %we will work in your root
fname = dir([root '/*.json']); %find all json files
fname = {fname.name}; 

for jsonnum = 1:length(fname) %go though your json files
    fid = fopen([root '/' fname{jsonnum}],'rt') ; %open json file
    disp(['Loaded JSON file!']);
    
    
    X = fread(fid) ; %read json file
    fclose(fid) ; %close handle
    X = char(X.') ; 
    disp(['Making changes to JSON file!']);
    
    %Making JSON readable for the Matlab decoder
    Y = strrep(X, 'ObjectId(', '') ;
    Y = strrep(Y, '),', ',') ;
    Y = strrep(Y, 'ISODate(', '') ;
    Y = strrep(Y, 'MongoDB shell version: 3.0.15', '') ;
    Y = strrep(Y, 'connecting to: omt', '') ;
    Y = regexprep(Y, '}[+\n\r]+{', '},{');
    
    Y = ['[' newline Y newline ']'];
    [~,nameonly,~] = fileparts([root '/' fname{jsonnum}]);
    %file name of corrected json file
    fname{jsonnum} = [nameonly '_corrected.json'];
    
    disp(['Saving new JSON file!']);
    
    %open new corrected json file
    fid2 = fopen([root '/' fname{jsonnum}],'wt') ;
    %save your text there
    fwrite(fid2,Y) ;
    fclose (fid2) ;
    
    clear X Y fid fid2
    
    disp(['Decoding JSON File!']);
    
    %decode your json file
    result = jsondecode(fileread([root '/' fname{jsonnum}]));
    %if there is actual data, run this
    if isstruct(result)
        %put data in correct format
        result2 = result.sessions;
        [result2.userName] = deal(result.x_id);
        [result2.task]     = deal(result2.versionName);
        result = result2;
        clear result2
        
        %=============
        % Find all Fractals 1&3 and completed files only
        indices = find((strcmp({result.task},'Fractals1&3') | strcmp({result.task},'Fractals_1&3') ...
            | strcmp({result.task},'Sanjay_Yoni')) & strcmp({result.state},'complete'));
        
        counter = 1;
        disp(['Creating mat files!']);
        
        %create all the variables of the very original fractals file format
        %and create a mat file that is identical to that format
        for i=indices
            time = [];
            block = [];
            trajectory_t = [];
            trajectory   = [];
            selectedId = [];
            targetId = [];
            identificationResponseTime = [];
            localisationResponseTime = [];
            endCoordinatesX = [];
            endCoordinatesY = [];
            foilId = [];
            targetOnTop = [];
            trial =[];
            allitemIds = [];
            allitemCoords = [];
            
            
            date = result(i).created;
            pid  = result(i).participantId;
            user = result(i).userName;
            condition = [];
            for j = 1:length(result(i).phases)
                for k = 1:length(result(i).phases(j).blocks)
                    block     = [block result(i).phases(j).blocks(k).trials.blockIndex];
                    trial     = [trial result(i).phases(j).blocks(k).trials.trialNumber];
                    condition = [condition {result(i).phases(j).blocks(k).trials.condition}];
                    
                    for m = 1:length(result(i).phases(j).blocks(k).trials)
                        time = [time {[result(i).phases(j).blocks(k).trials(m).screens(1).memoryStart ...
                            result(i).phases(j).blocks(k).trials(m).screens(1).memoryEnd ...
                            result(i).phases(j).blocks(k).trials(m).identification.identificationStart ...
                            result(i).phases(j).blocks(k).trials(m).identification.identificationResponseTime ...
                            result(i).phases(j).blocks(k).trials(m).identification.localisationResponseTime]}];
                        trajectory_t = [trajectory_t {cell2mat(cellfun(@(v)v(3),result(i).phases(j).blocks(k).trials(m).identification.track))'}];
                        trajectory   = [trajectory {complex(cell2mat(cellfun(@(v)v(1),result(i).phases(j).blocks(k).trials(m).identification.track)),cell2mat(cellfun(@(v)v(2),result(i).phases(j).blocks(k).trials(m).identification.track)))'}];
                        selectedId   = [selectedId result(i).phases(j).blocks(k).trials(m).identification.selectedId];
                        targetId     = [targetId result(i).phases(j).blocks(k).trials(m).identification.targetId];
                        identificationResponseTime = [identificationResponseTime result(i).phases(j).blocks(k).trials(m).identification.identificationResponseTime];
                        localisationResponseTime   = [localisationResponseTime result(i).phases(j).blocks(k).trials(m).identification.localisationResponseTime];
                        endCoordinatesX = [endCoordinatesX result(i).phases(j).blocks(k).trials(m).identification.endCoordinates(1)];
                        endCoordinatesY = [endCoordinatesY result(i).phases(j).blocks(k).trials(m).identification.endCoordinates(2)];
                        foilId = [foilId result(i).phases(j).blocks(k).trials(m).identification.nonTargetId];
                        targetOnTop = [targetOnTop result(i).phases(j).blocks(k).trials(m).identification.targetOnTop];
                        itemIds =[];
                        itemCoords = [];
                        for n = 1:length(result(i).phases(j).blocks(k).trials(m).screens)
                            itemIds =    [itemIds result(i).phases(j).blocks(k).trials(m).screens(n).items.index];
                            itemCoords = [itemCoords result(i).phases(j).blocks(k).trials(m).screens(n).items.coords];
                        end
                        itemCoords = itemCoords(:)';
                        
                        if size(itemIds,2)<size(allitemIds,2)
                            itemIds = [itemIds NaN(1,size(allitemIds,2)-size(itemIds,2))];
                            itemCoords = [itemCoords NaN(1,size(allitemCoords,2)-size(itemCoords,2))];
                            
                        elseif size(allitemIds,2)<size(itemIds,2)
                            allitemIds = [allitemIds NaN(size(allitemIds,1),size(itemIds,2)-size(allitemIds,2))];
                            allitemCoords = [allitemCoords NaN(size(allitemCoords,1),size(itemCoords,2)-size(allitemCoords,2))];
                        end
                        
                        
                        allitemIds    = [allitemIds; itemIds];
                        allitemCoords = [allitemCoords; itemCoords];
                    end
                    
                end
                
            end
            conditions = char(unique(condition));
            condition = cellfun(@(x)find(strcmp(x,unique(condition))),condition);
            data = [selectedId' targetId' identificationResponseTime' localisationResponseTime' endCoordinatesX' endCoordinatesY' allitemIds allitemCoords foilId' targetOnTop'];
            data_columns = {'selectedId','targetId','identificationResponseTime','localisationResponseTime','endCoordinatesX','endCoordinatesY'};
            
            %Item IDs
            for o = 1:size(allitemIds,2)
                data_columns = [data_columns {['itemId' num2str(o)]}];
            end
            
            for o = 1:size(allitemIds,2)
                data_columns = [data_columns {['location' num2str(o) 'X']} {['location' num2str(o) 'Y']}];
            end
            
            data_columns = char([data_columns {'foilId'} {'targetOnTop'}]);
            versionName = result(i).task;
            if strcmp(versionName,'Sanjay_Yoni') || strcmp(versionName,'Fractals_1&3') ...
                    || strcmp(versionName,'Sanjay_Yoni') || strcmp(versionName,'OMT_BioResources')
                versionName = 'Fractals1&3';
            end
            
            save([root '/' erase([versionName '_' result(i).userName '_' result(i).participantId '_' strrep(strrep(date,':','-'),'.','-') '.mat'],'/')],'block','condition','conditions','data','data_columns','date','pid','time','trajectory','trajectory_t','trial','user')
            clearvars -except result indices root counter i fname
            disp([num2str((counter/length(indices))*100) '% done']);
            counter = counter+1;
        end
    end
end

%go through all mat files
list = dir ([root '/*.mat']);
maxitems   = 4;                   %maximum amount of items that this dataset has
targetnum  = [1 1 1 1 1 1 1 1 1]; %target index for each condition
distractor = ...
    [99 99;              %distracter index for each condition, note 99 corresponds to the index of no distractor, in conditions where partcipants are presented with one fractal at encoding there ano no distractors presented
    99 99;
    2 3;
    2,3;
    99,99;
    99,99;
    99,99;
    2 3;
    2 3];


result       = [];
allfilenames = [];

%Loops through the list of all the datafiles.

k = 1;
m = 1;


for i=1:size(list,1)

    load([list(i).folder '/' list(i).name]) %loads each file
    result       = [result; data];          %adds each data to the variable result
    if any(condition==0)
        condition = condition+1;
    end
    cond(k,:)    = condition;               %adds the condition to an overall variable cond. So each line will be one subject
    [filename{1:length(condition)}] = deal(list(k).name);
    allfilenames = [allfilenames filename];
    k=k+1;
    excluded{m}=list(i).name;
    m=m+1;
    
end

Names = {list.name}';
Names = table(Names);

%==============================


%The variable result has all the data from all the participants:
%structure: each row will be one subject, each column will be a different parameter
allcond               = reshape(cond',1,size(cond,1)*size(cond,2));
alld                  = table2struct(array2table(result,'VariableNames',cellstr(data_columns)));
alld_diff             = reshape(alld,size(cond,2),size(cond,1))';
[alld(:).filename]    = deal(allfilenames{:});
allcond               = num2cell(reshape(cond',length(alld),1)');   %adds the condition number
[alld(:).condition]   = deal(allcond{:});
alltrialnum           = num2cell(reshape(repmat([1:size(cond,2)]',size(cond,1),1),length(alld),1))'; %adds to the trial number
[alld(:).trialnumber] = deal(alltrialnum{:});

%writetable(struct2table(alld), 'Results_preliminary.xlsx')

%% PREPARES DATA FOR ANALYSIS

%Reshapes the data from 'result'type into 'cond' type
%'result' has all the trials on it
%'cond' will have each line representing one participant and each column representing one trial
identification = [reshape([alld.selectedId]',size(cond,2),size(cond,1))' == reshape([alld.targetId]',size(cond,2),size(cond,1))'];
idnTime        = [reshape([alld.identificationResponseTime]',size(cond,2),size(cond,1))'];
locTime        = [reshape([alld.localisationResponseTime]',size(cond,2),size(cond,1))'];

failtrial      = (~identification); %where identified shape is not target shape
allcond        = cond; 
cartesianlocs  = [];

%extract locations, calculate distances to response
for i=1:maxitems
    location    = [[alld.(['location' num2str(i) 'X'])]' [alld.(['location' num2str(i) 'Y'])]'];
    response    = [[alld.endCoordinatesX]' [alld.endCoordinatesY]'];
    for j=1:size(location,1)
        dist(j)= pdist([location(j,:); response(j,:)]);
    end
    dist       = reshape(dist',size(cond,2),size(cond,1))';
    alldist(:,:,i)=dist;
    clear dist location response
end

for i=1:size(cond,1)
    for j=1:size(cond,2)
        targetdist(i,j)       = alldist(i,j,targetnum(cond(i,j))); %get the trial's target distance from response
        if ~any(distractor(cond(i,j),:)>maxitems)
            misbinditems(i,j) = min(alldist(i,j,distractor(cond(i,j),:)));
        else
            misbinditems(i,j) = NaN;
        end
        % John Grogan's Modelling =========== see documentation on this for
        % detail
        if cond(i,j)~=0
            cartesiantarg(i,j,:)   = [alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'X']);alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'Y'])];
            cartesianerr (i,j,:)   = [alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'X']);alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'Y'])] - [alld_diff(i,j).endCoordinatesX;alld_diff(i,j).endCoordinatesY];
            cartesianlocs_tmp = [];
            for k = 1:maxitems
                cartesianlocs_tmp = [cartesianlocs_tmp ; [alld_diff(i,j).(['location' num2str(k) 'X'])];[alld_diff(i,j).(['location' num2str(k) 'Y'])]];
            end
            cartesianlocs(i,j,:) = cartesianlocs_tmp;
            cartesiandist(i,j,:) = [sq(cartesianlocs(i,j,:)) - repmat(sq(cartesiantarg(i,j,:)),maxitems,1)];
        end
        %====================================
        response(i,j,:) = [alld_diff(i,j).endCoordinatesX alld_diff(i,j).endCoordinatesY];
    end
end

%Computation of Proportion of errors due to misbinding and guessing
%according to permutations. See "Donuts" paper for details
samplenum=5000;

for i=1:size(identification,1)
    for j=1:size(identification,2)
        misbindmatrix=[repmat(targetdist(i,j),samplenum,1) repmat(misbinditems(i,j),samplenum,1) datasample(rmmissing(misbinditems(i,:)),samplenum)'];
        [~,ii] = min(misbindmatrix,[],2);
        misbinding(i,j)  = sum(ii==2)/samplenum;
        guessrate(i,j)   = sum(ii==3)/samplenum;
        targetresp(i,j)  = sum(ii==1)/samplenum;
        imprecision(i,j) = min(misbindmatrix(1,1:2));
        clear misbindingmatrix ii
    end
end

misbinding_simu = NaN(size(cond,1),size(unique(allcond),1));
misbinding_MCMC = NaN(size(cond,1),size(unique(allcond),1));
guessrate_simu  = NaN(size(cond,1),size(unique(allcond),1));
guessrate_MCMC  = NaN(size(cond,1),size(unique(allcond),1));
precision_simu  = NaN(size(cond,1),size(unique(allcond),1));
precision_MCMC  = NaN(size(cond,1),size(unique(allcond),1));

for i = 1:15
    alldata_export{i} = NaN(size(cond,1),size(unique(cond),1));
end

distractor_nums = distractor;

%John Grogan's Modelling, see the Tutorial for Details
for i = 1:size(cond,1)
    k = 1;
    
    for condnum = [unique(allcond)]'
        currcond                              = condnum;
        distr_tmp                             = distractor_nums(currcond,:);
        distr_tmp(distr_tmp==99)              = targetnum(currcond);
        distractor_nums(currcond,:)           = distr_tmp;
        modeldata.items                       = rmmissing(sq(cartesianlocs(i,find(allcond(i,:)==condnum & ~failtrial(i,:)),:))');
        modeldata.errors                      = sq(cartesianerr(i,find(allcond(i,:)==condnum & ~failtrial(i,:)),:))';
        
        modeldata.distractors    = sq(cartesiandist(i,find(allcond(i,:)==condnum  & ~failtrial(i,:)),:))';
        if size(modeldata.distractors,1) == 1
            modeldata.distractors = modeldata.distractors';
        end
        modeldata.distractors    = modeldata.distractors(reshape([distractor_nums(currcond,:)*2-1; distractor_nums(currcond,:)*2], size(distractor_nums(currcond,:)*2-1,1),[]),:);
        modeldata.distractors(modeldata.distractors==0) = NaN;
        
        if sum(sum(~isnan(modeldata.distractors)))~=0 && size(modeldata.distractors,2)>4
            modeldata.distractors    = rmmissing(sq(modeldata.distractors),'MinNumMissing',size(modeldata.distractors,2));
            modeldata.dimensions     = [1366; 768];
            modeldata.minDists       = [88 29 88];
            model                    = SwapModel2D(); % the model to simulate
            modeldata.params         = [0.1, 0.1, 20]; % gamma, beta, SD
            tmp1 = MLE(modeldata, model);
            tmp2 = MemFit2D(modeldata, model, 'Verbosity', 0);
            misbinding_simu(i,k) = tmp1(1);
            misbinding_MCMC(i,k) = tmp2.posteriorMean(1);
            guessrate_simu(i,k)  = tmp1(2);
            guessrate_MCMC(i,k)  = tmp2.posteriorMean(2);
            precision_simu(i,k)  = tmp1(3);
            precision_MCMC(i,k)  = tmp2.posteriorMean(3);
            
           if ismember(condnum,cond)
               alldata_export{3}(i,condnum ) = tmp1(1);
               alldata_export{4}(i,condnum)  = tmp2.posteriorMean(1);
               alldata_export{6}(i,condnum)  = tmp1(2);
               alldata_export{7}(i,condnum)  = tmp2.posteriorMean(2);
               alldata_export{11}(i,condnum) = tmp1(3);
               alldata_export{12}(i,condnum) = tmp2.posteriorMean(3);
           end
           
           clear tmp1 tmp2 modeldata model

        end
        k = k+1;
    end
end

%put all data together
alldata_export{9}  = ones(size(alldata_export{3})) - (alldata_export{3} + alldata_export{6});
alldata_export{10} = ones(size(alldata_export{3})) - (alldata_export{4} + alldata_export{7});
alldata_export{1}  = groupMeans(targetdist+bool2nan(failtrial),2,cond);         %absolute error (pixels)
alldata_export{2}  = groupMeans(misbinding+bool2nan(failtrial),2,cond);         %proportion of errors due to misbinding
alldata_export{5}  = groupMeans(guessrate+bool2nan(failtrial),2,cond);          %proportion of errors due to guessing
alldata_export{8}  = groupMeans(targetresp+bool2nan(failtrial),2,cond);         %likelihood of choosing the target
alldata_export{13} = groupMeans(idnTime+bool2nan(failtrial),2,cond)/1000;       %identification time==time to identify the memoranda (seconds)
alldata_export{14} = groupMeans(locTime+bool2nan(failtrial),2,cond)/1000;       %localisation time== time to localise the memoranda (seconds)
alldata_export{15} = cellfun(@sum,groupMeans(identification,2,cond,@(x){x}))./cellfun('length',groupMeans(identification,2,cond,@(x){x})); %proportion of correctly identified trials
alldata_export{16} = groupMeans(imprecision+bool2nan(failtrial),2,cond);        %imprecision

%names of columns
DataNames = {'AbsoluteError_','Misbinding_Simple_','Misbinding_Simu_', ...
    'Misbinding_MCMC_','Guessing_Simple_','Guessing_Simu_', ...
    'Guessing_MCMC_','Target_Simple_','Target_Simu_', ...
    'Target_MCMC_','Precision_Simu_','Precision_MCMC_','IdentificationTime_','LocalisationTime_', ...
    'ProportionCorrect_','Imprecision_Simple_'};

%put together variable names
k = 1;
for i = 1:size(DataNames,2)
    for j = 1:size(conditions,1)
        allVariableNames{k} = [DataNames{i} strrep(strrep(strtrim(conditions(j,:)),' ','_'),',','') '_StarryNight'];
        k = k+1;
    end
end

%save all of this to csv file
EXPORT = horzcat(Names, array2table([alldata_export{:}],'VariableNames',allVariableNames));
writetable(EXPORT,[root '/YOUROMTCHART_' datestr(datetime('today')) '_' datestr(now,'HH_MM_SS_FFF') '.csv'])
