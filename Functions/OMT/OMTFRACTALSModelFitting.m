% OMTFRACTALSModelFitting
% Will load up the Fractals mat files created by OMTFRACTALS1AND3.m,
% and fit several models to them, and then do some figures and processing
% run in the same folder as OMTFRACTALS1AND3.m
% Requires matlib, ChangeMTB2DUnits.m from updated MemToolbox2D, and
% OMT2MemToolbox2D.m, and MakeRmanovaStats.m, and GetSubPlotShape.m

clear all;
%% set models

models = {  %StandardMixtureModel2D(),...
            SwapModel2D(),...
            SwapModelSepMisbind2D(),...
            WithRadialBias2D(SwapModel2D()),...
%             WithRadialBias2D(SwapModelSepMisbind2D),...
    };

nModels = length(models);
modelNames = {'Swap', 'Swap2', 'SwapBias'};%,'Swap2Bias'};
for i=1:nModels 
    nPars(1,i) = length(models{i}.paramNames); % number of parameters
end

%% what to split by?

% use {'isCorrect',1} only - for model comparison (fitting to all
% conditions at once)
fieldsToSplit = {   'isCorrect', 1; % only correct
                    'condition', 1:4; % main trials, no practice
                    }; % name of field, values to keep


minTrials = 20; % min number of trials needed for fitting - can lower if needed, but fewer trials = worse fitting

conds = {'1Short','1Long','3Short','3Long'}; % or change this to yours
% conds = {'correct'};
nConds = length(conds);

%% what units to use?

coordsUnits = 'mm'; % 'pixels', 'mm', or 'normalise'
% select your units here
fieldsToIgnore = {'condition','isPractice','isCorrect','pid','n','whichIsTestItem'}; % which fields to not convert units for
%% find files

resultFolder = './'; % where are your Fractals*.m files stored?
files = what(resultFolder);

fileString = 'Fractals'; % string for regexp

resFiles = files.mat(~cellfun(@isempty,regexp(files.mat, fileString)));

disp(resFiles); %%%%% the fittings/parameters below will be in the same order as these files

nFiles = length(resFiles); 


%% load files, split by condition

foilNumber = 4; % 4th item is the foil - will ignore
getPrevResps = 1; % will get previous response info

[data, dataByCond] = deal(cell(nFiles,1)); 
nTrials2 = NaN(nFiles,nConds); % to store number of trials per 
for i = 1:nFiles
    omtStruct = load(fullfile(resultFolder, resFiles{i})); % load the file
   
    data{i} = OMT2MemToolbox2D(omtStruct, foilNumber, getPrevResps); % format for memtoolbox2d
    data{i}.pid = omtStruct.pid; % store ID 

    % convert units?
    switch coordsUnits
        case 'pixels'
            % leave as is
        case 'normalise'
            data{i} = ChangeMTB2DUnits(data{i}, data{i}.dimensions/100,[],fieldsToIgnore); % express as % of screen dimensions
        case 'mm'
            data{i} = ChangeMTB2DUnits(data{i}, omtStruct.ppm,[],fieldsToIgnore); % express as % of screen dimensions
        otherwise
            error('coordsUnits must be pixels, noramlise or mm');
    end

    
    % split by conditions
    datasets = data{i};
    for iC = 1:size(fieldsToSplit, 1)
        datasets = SplitDataByField(datasets, fieldsToSplit{iC,1});
        check = cellfun(@(x) x.(fieldsToSplit{iC,1})(1), datasets); % get conditions
        toKeep = ismember(check, fieldsToSplit{iC,2}); % see which ones to keep
        if iC == size(fieldsToSplit,1)
            datasets = datasets(toKeep); % cell of 1 or many
        else
            datasets = datasets{toKeep}; % just one
        end
        
    end
    
    dataByCond{i} = datasets; % only keep those ones
    
    nTrials2(i,:) = cellfun(@(x) length(x.n),dataByCond{i});

end


%% fit models

% preallocate these variables to store the parameters - will work in parfor loop
clear fields fitPars negLogLike
for iC=1:nConds 
    for iM = 1:nModels
        fields{iM,iC} = ['model' char(iM+64) '_cond' num2str(iC)]; % get field names as placeholders for later
        
        fitPars.(fields{iM,iC}) = NaN(1,nPars(iM)); % pre fill with NaN in case of missing conditions
        negLogLike.(fields{iM,iC}) = NaN;
        
        me.(fields{iM,iC}) = {};
    end
end
fitPars = repmat(fitPars,nFiles,1); % repmat up to number of participants
negLogLike = repmat(negLogLike,nFiles,1);


% actually fit the models
for iPP = 1:nFiles %%%% switch to 'for' to 'parfor' if you have lots of files
    disp(iPP);
    for iM = 1:nModels
        for iC = 1:nConds
            if size(dataByCond{iPP}{iC}.n,2) >= minTrials
                [fitPars(iPP).(fields{iM,iC}),logLike] = MLE(dataByCond{iPP}{iC},models{iM}); % fit
                negLogLike(iPP).(fields{iM,iC}) = (-logLike); % logLike make negative            
            end    
        end
    end
end

%% extract stuff from structs

pars1 = NaN(nFiles,max(nPars),nModels,nConds); %preallocate
nll = NaN(nFiles,nModels,nConds);
for iC=1:nConds
    for iM=1:nModels
        pars1(:,1:nPars(iM),iM,iC) = nancat(1, fitPars.(fields{iM,iC}));
        nll(:,iM,iC) = nancat(1, negLogLike.(fields{iM,iC}));
    end
end

%% compare goodness of fit

% replace any Inf/-Inf
nll(isinf(abs(nll))) = NaN;

aic = 2 * (nll + nPars);
if isempty(nTrials2)
    bic = 2 * (nll) + nPars .* log(nTrials)';
else
    bic = 2 * (nll) + nPars .* log(permute(nTrials2,[1,3,2]));
end

disp(nanmean(bic,[1,3])); % mean BIC across pp + conds

[~,i] = min(nanmean(bic,[1 3]));
fprintf('the best fitting model is number %d: %s\n', i, models{i}.name);


%% reorder pars so it is [sd,a,g,b (or NaN),others]

nMaxPars = max(nPars);
pars = NaN(size(pars1));
pars(:,nMaxPars+1,:,:) = NaN;
for i=1:nModels
    parNames = models{i}.paramNames;
    
    sdInd = find(strcmpi(parNames,'sd')); % get sd params if there is an sd
    if sdInd
        sdPars = pars1(:,sdInd,i,:);
    else
        sdPars = [];
    end

    gInd = find(strcmpi(parNames,'g')); % find guessing param
    if gInd
        gPars = pars1(:,gInd,i,:);
    else
        gPars = [];
    end
    
    bInd = find(strcmpi(parNames,'b')); % find misbinding
    if bInd
        bPars = pars1(:,bInd,i,:); 
    else 
        bPars = []; 
    end
    
    if bInd % calculate target selection
        aPars = 1 - gPars - bPars;
    else
        aPars = 1-gPars;
    end
    
    % get any other params left - e.g. bias, misbinding2
    otherInds = find(~(strcmpi(parNames,'sd') | strcmpi(parNames,'g') | strcmpi(parNames, 'b')));
    if otherInds
        otherPars = pars1(:,otherInds,i,:);
    else
        otherPars = [];
    end
    
    % combine
    pars(:,1:nPars(i)+1,i,:) = [sdPars,aPars,gPars,bPars,otherPars];
    
end

% manually adjust some pars
pars(:,2,2,:) = 1 - sum(pars(:,3:5,2,:),2); % for sepmisbinding model, a=1-g-b1-b2


%% plot one model

for iModel = 1:3; % which to plot? - best fitting one

parNames = {'sd','a','g','b','mu'};% names of params in new order
nParams = max(nPars)+1; % number of params

factorLabels = {'1','3';'short','long'};
figure()
for j = iModel
    nP = nPars(iModel)+1;
    [r,c] = GetSubPlotShape(nP); % get rows and columns for subplotting
    for i = 1:nP
        par = permute(reshape(sq(pars(:,i,j,:)),nFiles,2,nConds/2),[1,2,3]); % permute as [1,3,2] to swap xaxis and lines
        if any(par,'all')
            subplot(r,c,i)
            h = errorBarPlot(par,'type','line',...
                'plotargs',{'LineWidth',2});
            set(gca,'XTick',1:4, 'XTickLabels',factorLabels(2,:));
            xlabel('condition')
            ylabel(parNames{i});
            
            box off
            xlim([0.75 2.25])
            
            
        end
    end
end
% SuperTitle(modelNames{j});
legend(factorLabels(1,:),'Location','Best')

end
[~,p] = ttest(sq(pars(:,5,2:3,:)),0,'tail','both');sq(p)

%% do anovas on each model

paramAn = cell(nModels,1);
paramStats = cell(nModels,1);
for j = 1:nModels
    for i = 1:nParams
        data = reshape(pars(:,i,j,:),nFiles,2,2);
        if any(data,'all')
            paramAn{j}{i} = rmanova(data, {'pp','delay','setSize'},'categorical',2:3);
        end
    end
    
    % make into table
    paramStats{j} = MakeRmanovaStats(paramAn{j},parNames(1:length(paramAn{j})));
    disp(paramStats{j});
end


