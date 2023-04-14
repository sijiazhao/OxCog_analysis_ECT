function data = OMT2MemToolbox2D(omtStruct, foilItemNumber, getPrevResps)
% Convert an OMT structure to a data structure for MemToolbox2D. Will have all the
% conditions/trials in one structure, so you can split these using
% SplitDataByField() or SplitDataByCondition(). 
% If omtStruct has a dimensions field it will use those dimensions,
% otherwise it will try to guess the dimensions based on item and response
% locations.
% Inputs:
%   omtStruct: matlab structure saved when processing the json file (see
%       OMTFRACTALS1AND3.m
%   foilItemNumber: item number of foil (if you want to remove it).
%       default=0 (i.e. none removed)
%   getPrevResps: 1 = get previous response, also store in biasCoords, and
%       distractors2. default is 0.
% 
% Outputs:
%   data = MemToolbox2D structure, containing fields: resps, items,
%   distLocs, targets, errors, conditions, n, dimensions, 
% 
% John Grogan, 2021

%% check inputs

if ~isstruct(omtStruct)
    error('omtStruct must be a structure');
end
if ~exist('foilItemNumber','var') || isempty(foilItemNumber)
    foilItemNumber = 0;
end
if ~exist('getPrevResps','var') || isempty(getPrevResps)
    getPrevResps = 0;
end
    

%% check fields exist

omtFields = fieldnames(omtStruct);
vitalFields = {'data','data_columns','condition'};

if ~all(ismember(vitalFields, omtFields))
    error('some of the vital fields (data, data_columns, condition) are missing');
end

dataCols = cellstr(omtStruct.data_columns);

nTr = size(omtStruct.data,1); % number of trials
%% get response

respCols = [find(strcmpi(dataCols, {'endCoordinatesX'})), find(strcmpi(dataCols, {'endCoordinatesY'}))];
data.resps = omtStruct.data(:,respCols)'; % repsonse coords [x;y]

%also store this in a cell for WithResponseSampling2D - (will not be affected by SplitDataByField)
data.respPdf = {data.resps};

%% get each item
itemNames = {'location%dX','location%dY'};

i = 1;
while true
    
    if i ~= foilItemNumber
        itemNamesI = {sprintf(itemNames{1}, i), sprintf(itemNames{2}, i)}; % make name
        if all(ismember(itemNamesI, dataCols)) % if exists
            [~, colInds] = ismember(itemNamesI, dataCols);
            data.items(i*2-1:i*2,:) = omtStruct.data(:, colInds)'; % item coords
        else % stop when both coords not found
            break
        end
    end
    
    i = i + 1;
     
end


%% remove bad items

% remove any that are all NaN
badItems = all(isnan(data.items), 2);
data.items(badItems,:) = [];

%% number of items per trial

data.n = sum(~isnan(data.items),1) ./ 2;

%% which is target

if isfield(omtStruct, 'whichIsTestItem')
    % use this
    
    data.whichIsTestItem = omtStruct.whichIsTestItem;
    n = 1:size(data.items,1);
    for i = 1:nTr % get this item for each trial
        inds = [data.whichItemIsTarget(i)*2-1; data.whichItemIsTarget(i)*2];
        data.targets(:,i) = data.items(inds,i);
        data.distLocs(:,i) = data.items(n~=inds,i);
    end
else
    data.whichIsTestItem = ones(1, nTr); % first item is test/target
    % assume 1st item is target
    data.targets = data.items(1:2,:);
    data.distLocs = data.items(3:end,:); % distractor lcoations
end

nDists = size(data.distLocs,1)/2;

%% errors
% resp - target
data.errors = data.resps - data.targets;

%% distractors
% dist - target
data.distractors = data.distLocs - repmat(data.targets, nDists,1) ;


%% other stuff

data.condition = omtStruct.condition; % condition

% practice conditions (or intro)
if isfield(omtStruct, 'conditions')
    pracInds = find(~cellfun(@isempty, regexp(cellstr(omtStruct.conditions), 'prac|intro')));
    
    data.isPractice = any(data.condition == pracInds,1); % is practice or intro
end


% isCorrectIdent
[isCorrInds, corrInds] = ismember({'selectedId', 'targetId'}, dataCols);
if all(isCorrInds)
    ids = omtStruct.data(:,corrInds);
    data.isCorrect = [ids(:,1) == ids(:,2)]'; % does selected match target? 
end

%% get previous response location?

if getPrevResps
    
    % previous trial's response location
    data.prevResp = [NaN(2,1), data.resps(:, 1:end-1)];

    data.biasCoords = data.prevResp; % for WithRadialBias2D
    data.distractors2 = data.prevResp; % for SwapModelSepMisbind2D
    data.distractors1 = data.distractors; % current trial distractors
end


%% dimensions

if isfield(omtStruct, 'dimensions')
    data.dimensions = omtStruct.dimensions;

else % will have to guess at dimensions
    
    if isfield(omtStruct, 'longDim') % get screen info
        scrDimMm = [omtStruct.longDim; omtStruct.longDim * .695]; % estimate y size from ratio of ipad (and of scrDimGuess ratio)
        scrDimPix = round(scrDimMm .* omtStruct.ppm); % in pixels
        scrDimMm = round(scrDimMm); % round it now
        
    else
        [scrDimMm, scrDimPix] = deal(NaN(2,1));
    end
    
    % get item + repsonse locations
    allItemLocs = reshape(data.items,2,[]); % [x;y]
    
    minMaxLocs = [min(allItemLocs,[],2), max(allItemLocs,[],2)]; % cols = x,y, rows = min, max. due to screen edge distances, the sum of these is close to the screen dimensions
    minMaxResps = round([min(data.resps,[],2), max(data.resps,[],2)]); % same for response
    
    scrDimGuess = max(minMaxLocs(:,1)) + minMaxLocs(:,2); % min = edge distance, max+min = top edge
    
    % pick screen dimensions
    if all(max(minMaxResps,[],2) <= scrDimGuess) % if all resps are within this screen, use that
        data.dimensions = scrDimGuess;
        disp('calculating screen dimensions based on item locations');
    
    elseif all(max(minMaxResps,[],2) <= scrDimPix) % if all within scrDimPix
        data.dimensions = scrDimPix;
        disp('calculating screen dimensions based on longDim + ppm');
    
    else % use max response
        data.dimensions = round(max(minMaxResps,[],2));
        disp('calculating screen dimensions based on max response locations');    
    end
    
    
end


