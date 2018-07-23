function [ output_args ] = struct_tree( StructIn, max_level )
%STRUCT_TREE Print out the hierarchy of a structure
%   Takes a structure as an input and prints out each field's name and - if
%   the structure is scalar - size and class. If that field is a structure
%   itself, it recursively prints out the subfields. Takes an optional
%   second argument, max_level that limits how deep it will recurse - e.g.
%   if this is set to 1, only the first level of fields will be printed.



E = JLLErrors;
narginchk(1,2)
if nargin < 2
    max_level = Inf;
elseif ~isnumeric(max_level) || ~isscalar(max_level)
    E.badinput('max_level must be a scalar number');
elseif max_level < 1
    E.badinput('max_level must be at least 1')
end

if ~isstruct(StructIn)
    E.badinput('Input must be a structure')
end

fprintf('%s\n',inputname(1));
check_for_fields(StructIn, 0, max_level);


end

function check_for_fields(StructIn, lvl, max_level)
lvl = lvl+1;
fns = fieldnames(StructIn);
spacer = repmat('    ',1,lvl);
for f=1:numel(fns)
    if numel(StructIn)==1
        size_str = mat2str(size(StructIn.(fns{f})));
        class_str = class(StructIn.(fns{f}));
    else
        size_str = '';
        class_str = '';
    end
    
    fprintf('%s%s - %s %s\n',spacer,fns{f},size_str,class_str);
    if isstruct(StructIn.(fns{f})) && lvl < max_level
        check_for_fields(StructIn.(fns{f}),lvl,max_level)
    end
end
end