function [ vec_out, std_err ] = avg_n_elements( V, n, varargin )
%AVG_N_ELEMENTS Averages together sets of n element in vector V
%   This function takes a vector V and will take every n elements,
%   calculate the average, median or mode of them, and return a vector of
%   length floor(length(V)/n). Requires a vector and number of elements,
%   also has several parameters.
%
%   'op' - the operation to perform. Defaults to 'mean', can be 'median',
%   'mode', 'nanmean', or 'nanmedian'.
%
%   'skip_end' - what to do with "extra" elements. If V is 23 elements long
%   and n is 5, by default the last three elements will just be ignored and
%   the returned vector will be 4 elements long. Set this to false to keep
%   those final elements.
%
%   'skip_start' - like skip_end, but for elements at the beginning. Does
%   nothing if start_ind is not modified also.
%
%   'start_ind' - allows the user to start the averaging later in the
%   vector. Useful if, for example, you want minute averages of 1 second
%   data, but the data starts in between minutes. Defaults to 1, i.e. does
%   not skip any data at the beginning.
%
%   Outputs the vector of means/medians/modes as well as a vector of
%   standard errors.
%
%   Josh Laughner <joshlaugh5@gmail.com> 12 June 2015

E = JLLErrors;
DEBUG_LEVEL = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

narginchk(2,Inf)

p = inputParser;
p.addParameter('op','mean');
p.addParameter('skip_start', true);
p.addParameter('skip_end', true);
p.addParameter('start_ind', 1);

p.parse(varargin{:});
pout = p.Results;

op_string = lower(pout.op);
skip_start = pout.skip_start;
skip_end = pout.skip_end;
start_ind = pout.start_ind;

if ~isvector(V) || ~isnumeric(V)
    E.badinput('V is expected to be a numeric vector')
elseif ~isscalar(n) || ~isnumeric(n)
    E.badinput('n is expected to be a numeric scalar')
elseif ~ismember(op_string,{'mean','nanmean','median','nanmedian','mode'})
    E.badinput('The parameter "op" must be mean, nanmean, median, nanmedian, or mode');
elseif ~isscalar(skip_start) || (~isnumeric(skip_start) && ~islogical(skip_start))
    E.badinput('The parameter "skip_start" is expected to be a boolean (numeric scalar)')
elseif ~isscalar(skip_end) || (~isnumeric(skip_end) && ~islogical(skip_end))
    E.badinput('The parameter "skip_end" is expected to be a boolean (numeric scalar)')
elseif ~isscalar(start_ind) || ~isnumeric(start_ind) || start_ind < 1 || start_ind > length(V)
    E.badinput('The parameter "start_ind" is expected to be a scalar with value between 1 and length(V)')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

V2 = V(start_ind:end);
last_ind = length(V) - mod(length(V2),n); % remove ending values such that we can make V2 into a matrix
V2 = V(start_ind:last_ind);

% Initialize the output vectors with enough room to hold the starting and
% ending orphans if requied.
num_avgs = floor(numel(V)/n) + ~skip_start + ~skip_end;
vec_out = zeros(1, num_avgs);
if iscolumn(V)
    vec_out = vec_out';
end
std_err = zeros(size(vec_out));

% Average (or median or mode) the bulk of the vector
Op = str2func(op_string);
V2 = reshape(V2, [n, numel(V2)/n]);
avgs = Op(V2,1);
% Handle whether the standard errors should ignore nans or not. This way if
% the user is expecting the presence of of NaNs to return one, it will be
% consistent.
if ~isempty(regexp(op_string,'nan','ONCE'))
    notnans = sum(~isnan(V2),1);
    serrs = nanstd(V2,1) ./ notnans;
    ignore_nans = true;
else
    serrs = std(V2,1) / n;
    ignore_nans = false;
end

% Put this in the output vectory in the right place
s = 1 + ~skip_start;
e = length(vec_out) - ~skip_end;
vec_out(s:e) = avgs;
std_err(s:e) = serrs;

% Handle the starting and ending bits if requested
if ~skip_start && start_ind > 1
    vec_out(1) = Op(V(1:start_ind));
    if ignore_nans
        std_err(1) = nanstd(V(1:start_ind));
    else
        std_err(1) = std(V(1:start_ind));
    end
end

if ~skip_end && last_ind < length(V)
    vec_out(end) = Op(V(last_ind+1:end));
    if ignore_nans
        std_err(end) = nanstd(V(last_ind+1:end));
    else
        std_err(end) = std(V(last_ind+1:end));
    end
end

