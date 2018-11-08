function [ cellout ] = sprintfmulti( formatstr, varargin )
%SPRINTFMULTI Expands FORMATSTR multiple times for array arguments
%   CELLOUT = SPRINTFMULTI( FORMATSTR, ... ) allows you to input a single
%   string with format markers as you would use with SPRINTF to insert
%   additional values, but this function will accept array arguments and
%   replicate FORMATSTR into a cell array of strings, CELLOUT, and
%   subsitute each value in the array in turn into each string.  This can
%   accept numeric or cell arrays so long as the contents of each cell in
%   the cell array is acceptable to SPRINTF.  Additional inputs that are
%   scalar numbers or individual strings will be inserted into every output
%   string.
%
%   Examples:
%
%   SPRINTFMULTI( '\\sigma = %.1f', [2, 3.6, 16.9] )
%   will return {'\sigma = 2.0', '\sigma = 3.6', '\sigma = 16.9'}.
%
%   SPRINTFMULTI( '%d: %s %s', 1:3, 'Alice', {'sits', 'stands', 'runs'})
%   will return {'1: Alice sits', '2: Alice stands', '3: Alice runs'}

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT CHECKING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;

if ~ischar(formatstr)
    E.badinput('FORMATSTR must be a string')
end

% Check that every value to insert into the string has the same number of
% elements or is scalar (and so should be put into every string). We have
% to check character arrays specially, because they have a length >= 1, but
% should be treated as a scalar.
n_el = nan(1, numel(varargin));
for i=1:numel(varargin)
    if ischar(varargin{i})
        n_el(i) = 1;
    else
        n_el(i) = numel(varargin{i});
    end
end
u_n_el = unique(n_el(n_el~=1));
if numel(u_n_el) > 1
    E.badinput('All additional array arguments to SPRINTFMULTI must have the same number of elements')
end

% If there are only scalar values, then the call to unique() will return an
% empty array, which causes repmat to issue a warning or error in versions
% of matlab after 2017.
if isempty(u_n_el)
    u_n_el = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%  MAIN FUNCTION  %%%%%
%%%%% (such as it is) %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

cellout = repmat({formatstr}, 1, u_n_el);
cellout = cprintf(cellout, varargin{:});

end

