function [ gobs, gnames ] = find_named_children( ax, varargin )
%FIND_NAMES_CHILDREN Finds named child graphics objects of the given axes
%   [ GOBS, GNAME ] = FIND_NAMED_CHILDREN( AX ) For the given axes, AX,
%   this will find all graphics objects that are children of AX and have a
%   non-empty DisplayName property. The handles to those objects will be
%   returned as the array GOBS and the DisplayNames as the cell array
%   GNAMES. This is particularly useful for recreating legends.
%
%   Parameter values:
%       'types' - a cell array of strings defining which types of graphics
%       objects should be returned. This works by filtering each child's
%       Type property against this cell array. Only objects whose Type is
%       contained in the cell array are returned.

E = JLLErrors;
p = inputParser;
p.addParameter('types', {'all'});
p.parse(varargin{:});
pout = p.Results;
child_types = pout.types;

if ischar(child_types)
    child_types = {child_types};
elseif ~iscellstr(child_types)
    E.badinput('CHILD_TYPES must be a cell array of strings');
end

ch = get(ax, 'children');
ch = flip(ch); % get almost always returns things backwards from the order they were created
named_bool = false(size(ch));
for a=1:numel(ch)
    if ~ismember('all', child_types) && ~ismember(ch(a).Type, child_types)
        continue
    elseif isempty(ch(a).DisplayName)
        continue
    end

    named_bool(a) = true;
end

gobs = ch(named_bool);
gnames = {ch(named_bool).DisplayName};

end

