function [ cell_out ] = glob( cell_in, pat )
%GLOB Applies regular expressions to a cell array.
%   In bash, a glob is what lets you type "ls *.m" to see all files ending
%   in .m - it's what expanding that * is called. This will do something
%   similar, filtering a cell array for strings that match the given
%   regular expression, pat.

E = JLLErrors;

%%%%% INPUT CHECKING %%%%%
if ~iscell(cell_in)
    E.badinput('First argument must be a cell array')
elseif any(~iscellcontents(cell_in,'ischar'))
    E.badinput('Cell array must contain only strings')
elseif ~ischar(pat)
    E.badinput('Second argument must be a string')
end

%%%%% GLOB %%%%%

xx = ~iscellcontents(regexp(cell_in,pat),'isempty');
cell_out = cell_in(xx);

end

