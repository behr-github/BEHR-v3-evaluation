function [ cellout ] = cprintf( formatcell, varargin )
%CELLOUT = CPRINTF( FORMATCELL, ... ) Format each string in FORMATCELL
%   Essentially a version of SPRINTF that accepts a cell array of format
%   strings and will insert all the following arguments into each string
%   using sprintf.
%
%   If there are format arguments whose values you want to change for each
%   cell in CELLOUT, pass their values as either arrays or cell arrays. For
%   example:
%
%   cprintf({'%d: %s %s is unknown', '%d: %s %s is scary'},[1 2],'John',{'Doe','Reese'})
%
%   will return:
%
%   CELLOUT = {'1: John Doe', '2: John Snow'}

E = JLLErrors;
if ~iscellstr(formatcell)
    E.badinput('FORMATCELL must be a cell array of strings');
end

for i=1:numel(formatcell)
    % Construct a temporary cell array that holds the corresponding
    % elements of the arguments
    argcell = cell(size(varargin));
    for j=1:numel(varargin)
        if isscalar(varargin{j}) || ischar(varargin{j})
            argcell{j} = varargin{j};
        elseif iscell(varargin{j})
            argcell{j} = varargin{j}{i};
        else
            argcell{j} = varargin{j}(i);
        end
    end
    try
        formatcell{i} = sprintf(formatcell{i}, argcell{:});
    catch err
        % Print out a slightly more helpful error message so that the user
        % knows which cell to look at.
        error(err.identifier, 'Problem with cell %d:\n%s',i,err.message);
    end
end

cellout = formatcell;

end

