function varargout = parsepvpairs(names, defaults, varargin)
%PARSEPVPAIRS Validate parameter name/value pairs and throw errors if necessary.
%   Given the cell array of valid parameter names, a corresponding cell array of
%   parameter default values, and a variable length list of parameter name/value
%   pairs, validate the specified name/value pairs and assign values to output
%   parameters.
%
%   [P1, P2, ...] = parsepvpairs(Names, Defaults, 'Name1', Value1, 'Name2', Value2, ...)
%
%   Inputs:
%      Names - Cell array of valid parameter names.
%
%   Defaults - Cell array of default values for the parameters named in Names.
%
%      Name# - Character strings of parameter names to be validated and
%              assigned the corresponding value that immediately follows each
%              in the input argument list. Parameter name validation is
%              case-insensitive and partial string matches are allowed provided
%              no ambiguities exist.
%
%     Value# - The values assigned to the corresponding parameter that
%              immediately precede each in the input argument list.
%
%   Outputs:
%         P# - Parameters assigned the parameter values Value1, Value2, ...
%              in the same order as the names listed in Names. Parameters
%              corresponding to entries in Names that are not specified in the
%              name/value pairs are set to the corresponding value listed in
%              Defaults.

%   Copyright 1995-2007 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $   $Date: 2008/12/21 01:51:18 $

%
% Short-circuit the input argument checking for performance purposes. Under the
% following specific circumstances, users may by-pass P-V pair validation when:
%
% (1) Values are assigned to all recognized parameters,
% (2) Parameters are specified in exactly the same order as in the input NAMES,
% (3) All parameters are completely specified (i.e., exact matches and no partial 
%     names allowed).
%

if isequal(varargin(1:2:end-1), names)
   varargout = varargin(2:2:end);
   return
end 

% Initialize some variables.
nInputs   = length(varargin);  % # of input arguments
varargout = defaults;

% Ensure parameter/value pairs.
if mod(nInputs, 2) ~= 0
   error('finance:parsepvpairs:incorrectNumberOfInputs', ...
      'Input parameters must be in name/value pair format.');

else
   names = lower(names);

   % Process p/v pairs.
   for j = 1:2:nInputs
      pName = varargin{j};

      if ~ischar(pName)
         error('finance:parsepvpairs:nonTextString', ...
            'Parameter names must be character strings.');
      end

      i = strmatch(lower(pName), names);

      if isempty(i)
         error('finance:parsepvpairs:invalidParameter', ...
            'Invalid parameter name:  %s.', pName);

      elseif length(i) > 1
         % If ambiguities exist, check for exact match to narrow search.
         i = strmatch(lower(pName), names, 'exact');
         if length(i) == 1
            varargout{i} = varargin{j+1};

         else
            error('finance:parsepvpairs:ambiguousParameter', ...
               'Ambiguous parameter name:  %s.', pName);
         end

      else
         varargout{i} = varargin{j+1};
      end
   end
end


% [EOF]