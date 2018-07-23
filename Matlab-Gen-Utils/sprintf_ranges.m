function [ str ] = sprintf_ranges(values, varargin)
%SPRINTF_RANGES Print a string representing ranges of values.
%   STR = SPRINTF_RANGES( VALUES ) Returns a character array, STR, that
%   lists ranges of values with a difference of 1 between each in VALUES
%   separated by a dash and separate ranges separated by an underscore. For
%   example:
%
%       sprintf_ranges(1:10) = '1-10'
%       sprintf_ranges([1, 2, 3, 4, 5, 6]) = '1-6'
%       sprintf_ranges([1, 2, 3, 10, 18, 19, 20]) = '1-3_10_18-20'
%
%   In the first two examples, all the values 1 through 10 or 1 through 6
%   are separated by 1, so SPRINTF_RANGES assumes that they are part of a
%   contiguous range and prints them as such. In the third example, there
%   are two contiguous ranges, 1 through 3 and 18 through 20, and one
%   "orphan" value, 10. Each group (range or orphan) is printed, separated
%   from the others by an underscore.
%
%   STR = SPRINTF_RANGES( VALUES, INCREMENT ) The value for INCREMENT
%   alters the default assumed increment in continuous ranges. So, for
%   example,
%
%       sprintf_ranges([2, 4, 6, 8]) = '2_4_6_8'
%       sprintf_ranges([2, 4, 6, 8], 2) = '2-8
%
%   In the first example, since each value is separated by 2 from its
%   neighbors, SPRINTF_RANGES assumes that they are all "orphans". In the
%   second example, including INCREMENT = 2 tells it that contigous ranges
%   should be separated by 2, and so it collapses the input into a single
%   range.
%
%   Additional parameters:
%
%       'range_sep': the character(s) placed between the two ends of a
%       range. Default is '-'.
%
%       'value_sep': the character(s) placed between separate groups
%       (ranges or orphans). Default is '_'.
%
%       'rel_tol': the relative tolerance for determining whether a
%       difference is equal to the increment. Simply checking difference ==
%       increment is prone to floating point errors, so instead this
%       function checks:
%
%           abs( difference - increment ) < rel_tol * increment
%
%       The default value for rel_tol is 1e-4, i.e. differences need only
%       be within 0.01% of the increment to be considered equal.

E = JLLErrors;

p = inputParser;
p.addOptional('increment',1,@isnumeric);
p.addParameter('range_sep','-');
p.addParameter('value_sep','_');
p.addParameter('rel_tol',1e-4);

p.parse(varargin{:});
pout = p.Results;

increment = pout.increment;
range_sep = pout.range_sep;
value_sep = pout.value_sep;

if ~isnumeric(values)
    E.badinput('VALUES must be numeric');
end

if ~isnumeric(increment) || ~isscalar(increment)
    E.badinput('INCREMENT must be a scalar number');
end
if ~ischar(range_sep)
    E.badinput('The parameter "range_sep" must be a character array');
end
if ~ischar(value_sep)
    E.badinput('The parameter "value_sep" must be a character array');
end
if ~isnumeric(pout.rel_tol) || ~isscalar(pout.rel_tol) || pout.rel_tol < 0
    E.badinput('The parameter "rel_tol" must be a positive scalar number');
end

tolerance = pout.rel_tol * increment;

% Main Function %

% Corner case: only one value passed. diff(values(:)) will then be an empty
% array.
if numel(values) == 1
    str = num2str(values);
    return
end

value_diffs = diff(values(:));
if all(mod(values,1) == 0)
    single_format = '%d';
    range_format = '%d%s%d';
else
    single_format = '%g';
    range_format = '%g%s%g';
end


i_diff = 1;
last_range_start = 1;

range_strs = {};

% get_next_group() will increment i_diff by variable amounts each time
% through the loop.
while i_diff <= numel(value_diffs)
    % This will get the last index of a group of values that are all
    % separated by the given increment.
    range_end = get_next_group();
    
    % If there are multiple adjacent values all separated by the requested
    % increment, group them into a range separated by the range separator.
    % Otherwise, just put the lone value by itself.
    if range_end == last_range_start
        range_strs{end+1} = sprintf(single_format, values(range_end)); %#ok<AGROW>
    else
        range_strs{end+1} = sprintf(range_format, values(last_range_start), range_sep, values(range_end)); %#ok<AGROW>
    end
    last_range_start = i_diff + 1;
    i_diff = i_diff + 1;
end

% If the last value given is part of a range separated by the requested
% increment, then the above loop will include it as the end of the final
% range. But, if it should be by itself, then it won't be captured, because
% there's one less difference than values. This block catches those orphan
% trailing values.
if abs(value_diffs(end) - increment) >= tolerance
    range_strs{end+1} = sprintf(single_format, values(i_diff));
end

str = strjoin(range_strs, value_sep);

    function incr_end = get_next_group()
        % This will find the next block of values that are each separated
        % by the given increment
        while i_diff <= numel(value_diffs) && abs(value_diffs(i_diff) - increment) < tolerance
            i_diff = i_diff + 1;
        end
        incr_end = i_diff;
    end


end
