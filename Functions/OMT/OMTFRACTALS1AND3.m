function OMTFRACTALS1AND3
%==================================================================
%A script put together by Younes Adam Tabi & Maria Raquel Maio
%==================================================================
% requires ChangeMTB2DUnits.m (in new release of MemToolbox2D) and matlib
% first half of script will convert json to mat files and save them. second
% half will do the analysis and save metrics/parameters
% 
% There are options to in the second half to deal with the screen size and
% model fitting (search 'ScreenSize')


%%

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
       try % this is only need for data from the server?
        result2 = result.sessions;
        [result2.userName] = deal(result.x_id);
        [result2.task]     = deal(result2.versionName);
        result = result2;
        clear result2
       end 
        %=============
        % Find all Fractals 1&3 and completed files only
        indices = find((strcmp({result.task},'Fractals1&3') | strcmp({result.task},'Fractals_1&3') ...
            | strcmp({result.task},'Sanjay_Yoni')) & strcmp({result.state},'complete'));
        
        counter = 1;
        disp(['Creating mat files!']);
        
        %create all the variables of the very original fractals file format
        %and create a mat file that is identical to that format
        for i=indices % each file
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
            for j = 1:length(result(i).phases) % phase
                for k = 1:length(result(i).phases(j).blocks) % block
                    block     = [block result(i).phases(j).blocks(k).trials.blockIndex];
                    trial     = [trial result(i).phases(j).blocks(k).trials.trialNumber];
                    condition = [condition {result(i).phases(j).blocks(k).trials.condition}];
                    
                    for m = 1:length(result(i).phases(j).blocks(k).trials) % trial
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
            for o = 1:size(allitemIds,2) % item (incl foil)
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
            
            % insert screen dimension info
            longDim = result(i).meta.calibration.longDimension; % length of longest screen dimension
            ppm = result(i).meta.calibration.ppm; % number of pixels per mm
            scaleFactor = result(i).meta.calibration.scaleFactor; % scaling factor of items
            rngSeed = result(i).meta.rngSeed; % random seed for stim locations
            
            save([root '/' erase([versionName '_' result(i).userName '_' result(i).participantId '_' strrep(strrep(date,':','-'),'.','-') '.mat'],'/')],'block','condition','conditions','data','data_columns','date','pid','time','trajectory','trajectory_t','trial','user','longDim','ppm','scaleFactor','rngSeed')
            clearvars -except result indices root counter i fname
            disp([num2str((counter/length(indices))*100) '% done']);
            counter = counter+1;
        end
    end
end

%go through all mat files
list = dir ([root '/Fractals*.mat']);
maxitems   = 4;                   %maximum amount of items that this dataset has (incl foil)
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

scr = struct();
for i=1:size(list,1) % each file
    
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
    
    scr(i).longDim = longDim; % long dimension (mm)
    scr(i).ppm = ppm; % pixels per mm
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
for i=1:maxitems % each item
    location    = [[alld.(['location' num2str(i) 'X'])]' [alld.(['location' num2str(i) 'Y'])]'];
    response    = [[alld.endCoordinatesX]' [alld.endCoordinatesY]'];
    for j=1:size(location,1)
        dist(j)= pdist([location(j,:); response(j,:)]); % euclidean distance of item - response
    end
    dist       = reshape(dist',size(cond,2),size(cond,1))';
    alldist(:,:,i)=dist;
    clear dist location response
end

for i=1:size(cond,1) % participant
    for j=1:size(cond,2) % trial
        targetdist(i,j)       = alldist(i,j,targetnum(cond(i,j))); %get the trial's target distance from response
        if ~any(distractor(cond(i,j),:)>maxitems) % if there are nonTargers
            misbinditems(i,j) = min(alldist(i,j,distractor(cond(i,j),:))); % smallest distance to nonTarget
        else
            misbinditems(i,j) = NaN;
        end
        % John Grogan's Modelling =========== see documentation on this for
        % detail
        if cond(i,j)~=0
            cartesiantarg(i,j,:)   = [alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'X']);alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'Y'])]; % coords of target [pp trial x/y]
            cartesianerr (i,j,:)   = [alld_diff(i,j).endCoordinatesX;alld_diff(i,j).endCoordinatesY] - [alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'X']);alld_diff(i,j).(['location' num2str(targetnum(cond(i,j))) 'Y'])]; % response - target
            cartesianlocs_tmp = [];
            for k = 1:maxitems % all items
                cartesianlocs_tmp = [cartesianlocs_tmp ; [alld_diff(i,j).(['location' num2str(k) 'X'])];[alld_diff(i,j).(['location' num2str(k) 'Y'])]];
            end
            cartesianlocs(i,j,:) = cartesianlocs_tmp; % all locs [pp trial [item1x, item1y, item2x...]]
            cartesiandist(i,j,:) = [sq(cartesianlocs(i,j,:)) - repmat(sq(cartesiantarg(i,j,:)),maxitems,1)]; % item - target
        end
        %====================================
        response(i,j,:) = [alld_diff(i,j).endCoordinatesX alld_diff(i,j).endCoordinatesY]; % response location
    end
end

%Computation of Proportion of errors due to misbinding and guessing
%according to permutations. See "Donuts" paper for details
samplenum=5000;

for i=1:size(identification,1) % participant
    for j=1:size(identification,2) % trial
        misbindmatrix=[repmat(targetdist(i,j),samplenum,1) repmat(misbinditems(i,j),samplenum,1) datasample(rmmissing(misbinditems(i,:)),samplenum)']; % [samplenum x 3] matrix of [targetDistance, nearest nonTarget distance, random trials' misbinding distances]
        [~,ii] = min(misbindmatrix,[],2); % which column is smallest
        misbinding(i,j)  = sum(ii==2)/samplenum; % mean rate of nonTargDistance being smallest
        guessrate(i,j)   = sum(ii==3)/samplenum; % mean rate of random trial nonTargDistance
        targetresp(i,j)  = sum(ii==1)/samplenum; % mean rate of target distance being smallest
        imprecision(i,j) = min(misbindmatrix(1,1:2)); % nearest neighbour distance
        clear misbindingmatrix ii
    end
end

misbinding_simu = NaN(size(cond,1),size(unique(allcond),1));% preallocate these for storing modelling params
misbinding_MCMC = NaN(size(cond,1),size(unique(allcond),1));
guessrate_simu  = NaN(size(cond,1),size(unique(allcond),1));
guessrate_MCMC  = NaN(size(cond,1),size(unique(allcond),1));
precision_simu  = NaN(size(cond,1),size(unique(allcond),1));
precision_MCMC  = NaN(size(cond,1),size(unique(allcond),1));

for i = 1:17 % preallocate this for also storing params
    alldata_export{i} = NaN(size(cond,1),size(unique(cond),1));
end

distractor_nums = distractor;

%% set ScreenSize and model fitting options here

%%%%%% if all your people have the same screen size, set this to 1
allSameScreenSize = 0 % do they all have the same screen size? [1 or 0]
if allSameScreenSize %%%%% and enter the screen size here
    screenSize = [2048; 1536]; % [x;y] size in pixels - this is the ipad air 2 (see FAQ on OMT)
    coordsUnits = 'pixels'; % can still use 'normalised' or 'mm'
    warning('all participants used the SAME size screen, so using size [%d;%d] in %s\n', screenSize(1), screenSize(2), coordsUnits);
else
    % if they have different screen sizes, this will affect the scaling of
    % the imprecision parameter (and potentially the others) meaning
    % between subject effects may be due to differences in screens used. so
    % it is recommended to convert all coordinates to % of screen or
    % millimetres
    coordsUnits = 'normalise'; % or 'mm'. pixels not recommended
    warning('participants used DIFFERENT size screens, so will be fitting with units: %s\n', coordsUnits);
end

% do you want to run MCMC fitting as well as MLE? (takes much longer)
doMCMC = 0; % 1 = use MemFit2D to to MCMC fitting (currently only the posterior mean is saved)


%%%%% minimum number of trials needed for each condition for fitting?
minTrials = 20; % 20 is best. 10 is ok. 4 is absolute minimum needed (will be very noisy)

%% John Grogan's Modelling, see the Tutorial for Details
for i = 1:size(cond,1) % pp
    k = 1; %condnum, but only for the 3-item main trials
    
    if allSameScreenSize % if they all use the same
        
        dimensions = screenSize; % those are the dimensions
        
    else % if using different screen sizes, will attempt to calculate the 
         % dimensions, which can be used to normalise later
        
        % set screen dimensions based on longDim, ppm
        scrDimMm = [scr(i).longDim; scr(i).longDim * .695]; % estimate y size from ratio of ipad (and of scrDimGuess ratio)
        scrDimPix = round(scrDimMm .* scr(i).ppm); % in pixels
        scrDimMm = round(scrDimMm); % round it now

        % the items were a certain distance from the screen edges, so
        % estimate that
        allItemLocs = reshape(sq(cartesianlocs(i,:,:))',2,[]); % get all item locations for this person, all conds.
        minMaxLocs = [min(allItemLocs,[],2), max(allItemLocs,[],2)]; % rows = x,y, cols = min, max. due to screen edge distances, the sum of these is close to the screen dimensions
        minMaxResps = [min(sq(response(i,:,:)))', max(sq(response(i,:,:)))']; % same for response

        % estimate the dimensions
        scrDimGuess = max(minMaxLocs(:,1)) + minMaxLocs(:,2); % min = edge distance, max+min = top edge
        
        % pick screen dimensions - different options
        if all(max(minMaxResps,[],2) <= scrDimGuess) % if all resps are within this screen, use that
            dimensions = scrDimGuess;
            disp('calculating screen dimensions based on item locations');
        elseif all(max(minMaxResps,[],2) <= scrDimPix) % if all within scrDimPix
            dimensions = scrDimPix;
            disp('calculating screen dimensions based on longDim + ppm');
        else % use max response
            dimensions = round(max(minMaxResps,[],2));
            disp('calculating screen dimensions based on max response locations');
        end
    end
    
    for condnum = [unique(allcond)]' %conds
        modeldata.dimensions = dimensions; % use the ones set above
        
        currcond                              = condnum;
        distr_tmp                             = distractor_nums(currcond,:); % get number of distractors
        distr_tmp(distr_tmp==99)              = targetnum(currcond); % replace 99 with 1
        distractor_nums(currcond,:)           = distr_tmp; % store this into distractors_nums
        modeldata.items                       = rmmissing(sq(cartesianlocs(i,find(allcond(i,:)==condnum & ~failtrial(i,:)),:))'); % only trials correctly identified, remove NaN distractors
        modeldata.errors                      = sq(cartesianerr(i,find(allcond(i,:)==condnum & ~failtrial(i,:)),:))'; % resp-target
        modeldata.distractors    = sq(cartesiandist(i,find(allcond(i,:)==condnum  & ~failtrial(i,:)),:))'; % keep distractors in here. item - target
        
        if size(modeldata.distractors,1) == 1 % if there was only one item, the sq() above makes this into 1 row, so rotate it
            modeldata.distractors = modeldata.distractors';
        end
        % i think this bit removes foils that were put into distractors -
        modeldata.distractors    = modeldata.distractors(reshape([distractor_nums(currcond,:)*2-1; distractor_nums(currcond,:)*2], size(distractor_nums(currcond,:)*2-1,1),[]),:); % keep only distractors?
        modeldata.distractors(modeldata.distractors==0) = NaN; % any where targ-item == 0 is target, so set to NaN
        
        if sum(~isnan(modeldata.errors),'all')>minTrials*2%sum(sum(~isnan(modeldata.distractors)))~=0 && size(modeldata.distractors,2)>4 % if there are more than 2 distractors. Change this - if nTrials > minNeeded
            % just leave as NaN
            %modeldata.distractors    = rmmissing(sq(modeldata.distractors),'MinNumMissing',size(modeldata.distractors,2)); % remove missing?
            
            % check units here, normalise or convert as needed
            switch coordsUnits
                case 'pixels'
                    % leave as is
                case 'normalise'
                    modeldata = ChangeMTB2DUnits(modeldata, modeldata.dimensions/100); % express as % of screen dimensions
                case 'mm'
                    modeldata = ChangeMTB2DUnits(modeldata, scr(i).ppm); % express as % of screen dimensions
                otherwise
                    error('coordsUnits must be pixels, noramlise or mm');
            end
            
            model = SwapModel2D(); % the model to simulate
            % params are [guess, misbind, imprecision]
                        
            [mlePars, mleLogLike] = MLE(modeldata, model); % maximum likelihood estimation of model. need to also get loglike
            % store params here (seems to be unused?)
            misbinding_simu(i,k) = mlePars(2);
            guessrate_simu(i,k)  = mlePars(1);
            precision_simu(i,k)  = mlePars(3);
            
            if doMCMC % do MCMC fitting - takes much longer per person
                % only really useful if you use the posteriors (i.e. not just
                % take the mean like is done here)
                tmp2 = MemFit2D(modeldata, model1); % Markov-Chain Monte Carlo fitting
                misbinding_MCMC(i,k) = tmp2.posteriorMean(2); % mean of posterior
                guessrate_MCMC(i,k)  = tmp2.posteriorMean(1);
                precision_MCMC(i,k)  = tmp2.posteriorMean(3);
            end
            
            % also store params here - exported below
            if ismember(condnum,cond) % if it is in conds (has to be?), store here for printing
                alldata_export{3}(i,condnum ) = mlePars(2); % misbinding
                alldata_export{6}(i,condnum)  = mlePars(1); % guessing
                alldata_export{11}(i,condnum) = mlePars(3);
                %store BIC
                alldata_export{17}(i, condnum) = -2*mleLogLike + log(size(modeldata.errors,2))*length(model.paramNames);
                
                if doMCMC % store these - NaN if not used
                    alldata_export{4}(i,condnum)  = tmp2.posteriorMean(2);
                    alldata_export{7}(i,condnum)  = tmp2.posteriorMean(1);
                    alldata_export{12}(i,condnum) = tmp2.posteriorMean(3);
                end
            end
            
            clear tmp1 tmp2 modeldata model
            
        end
        k = k+1;
    end
end

%put all data together - cell array of [pp x cond] means
alldata_export{9}  = ones(size(alldata_export{3})) - (alldata_export{3} + alldata_export{6}); % target rate = 1 - guessing - misbinding
alldata_export{10} = ones(size(alldata_export{3})) - (alldata_export{4} + alldata_export{7}); % same for MCMC
alldata_export{1}  = groupMeans(targetdist+bool2nan(failtrial),2,cond);         %absolute error (pixels)
alldata_export{2}  = groupMeans(misbinding+bool2nan(failtrial),2,cond);         %proportion of errors due to misbinding
alldata_export{5}  = groupMeans(guessrate+bool2nan(failtrial),2,cond);          %proportion of errors due to guessing
alldata_export{8}  = groupMeans(targetresp+bool2nan(failtrial),2,cond);         %likelihood of choosing the target
alldata_export{13} = groupMeans(idnTime+bool2nan(failtrial),2,cond)/1000;       %identification time==time to identify the memoranda (seconds)
alldata_export{14} = groupMeans(locTime+bool2nan(failtrial),2,cond)/1000;       %localisation time== time to localise the memoranda (seconds)
alldata_export{15} = cellfun(@sum,groupMeans(identification,2,cond,@(x){x}))./cellfun('length',groupMeans(identification,2,cond,@(x){x})); %proportion of correctly identified trials
alldata_export{16} = groupMeans(imprecision+bool2nan(failtrial),2,cond);        %imprecision (nearest neighbour distance)


%names of columns
DataNames = {'AbsoluteError_','Misbinding_Simple_','Misbinding_Simu_', ...
    'Misbinding_MCMC_','Guessing_Simple_','Guessing_Simu_', ...
    'Guessing_MCMC_','Target_Simple_','Target_Simu_', ...
    'Target_MCMC_','Precision_Simu_','Precision_MCMC_','IdentificationTime_','LocalisationTime_', ...
    'ProportionCorrect_','Imprecision_Simple_','BIC_Simple_'};

%put together variable names
k = 1;
for i = 1:size(DataNames,2) % names
    for j = 1:size(conditions,1) % each cond
        allVariableNames{k} = [DataNames{i} strrep(strrep(strtrim(conditions(j,:)),' ','_'),',','') '_StarryNight'];
        k = k+1;
    end
end

%save all of this to csv file
EXPORT = horzcat(Names, array2table([alldata_export{:}],'VariableNames',allVariableNames));
writetable(EXPORT,[root '/YOUROMTCHART_' datestr(datetime('today')) '_' datestr(now,'HH_MM_SS_FFF') '.csv'])
