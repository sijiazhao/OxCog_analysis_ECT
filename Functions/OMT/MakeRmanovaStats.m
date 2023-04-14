function statsTab = MakeRmanovaStats(stats, varNames)
% function statsTab = MakeRmanovaStats(stats, varNames)
% convert a cell array of rmanova statistics into a table
% 
% Input: stats = cell vector of rmanova outputs
%        varNames = optional cell array of names for table columns
% 
% Output: statsTab = table with row per term (excl intercept), column per
%                   variable
% 
n = length(stats); % number of cells

for i = 1:n
    
    pVals(:,i) = stats{i}.pValue(2:end); % get pvalues
        
end

% put into table, with rowNames
statsTab = array2table(pVals, 'RowNames', stats{1}.Term(2:end));

% varNames as columns
if exist('varNames','var') && ~isempty(varNames)
    statsTab.Properties.VariableNames = varNames;
end
