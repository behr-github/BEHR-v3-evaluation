function [ s_string ] = struct2string( S, varargin )
%STRUCT2STRING Returns a string representation of a structure.
%   S_STRING = STRUCT2STRING( S ) Returns a string representation of
%   structure S that gives the name of each field and its value or its type
%   and size (e.g. [1x2 cell]).
%
%   Parameters:
%
%       'sep' - separator character(s) that is/are printed between each
%       field. Default is a space.
%
%       'array_max' - maximum number of elements an array can have and be
%       printed explicitly before it is just printed as a type and size.
%       Default is 9.

E = JLLErrors;
p = inputParser;


p.addParameter('sep', ', ');
p.addParameter('array_max', 9);
p.parse(varargin{:});
pout = p.Results;

field_sep = pout.sep;
array_max_elements = pout.array_max;

if ~isstruct(S)
    E.badinput('S must be a structure');
elseif ~isscalar(S)
    E.notimplemented('%s','STRUCT2STRING is not implemented for non-scalar structures yet');
end

if ~ischar(field_sep)
    E.badinput('The parameter "field_sep" must be given a string');
end

if ~isnumeric(array_max_elements) || ~isscalar(array_max_elements) || array_max_elements < 0
    E.badinput('The parameter "array_max_elements" must be a scalar positive number');
end

fns = fieldnames(S);
val_strings = cell(size(fns));

for f=1:numel(fns)
    if ischar(S.(fns{f}))
        val = S.(fns{f});
    elseif isnumeric(S.(fns{f}))
        if isscalar(S.(fns{f}));
            val = sprintf('%g', S.(fns{f}));
        elseif numel(S.(fns{f})) < array_max_elements
            val = mat2str(S.(fns{f}));
        else
            val = format_typestr(S.(fns{f}));
        end
    else
        val = format_typestr(S.(fns{f}));
    end
    val_strings{f} = sprintf('%s: %s', fns{f}, val);
end

s_string = strjoin(val_strings, field_sep);

end

function tstr = format_typestr(val)
sz = num2cell(size(val));
sz_str = cellfun(@num2str, sz, 'uniformoutput', false);
sz_str = strjoin(sz_str, 'x');
tstr = sprintf('[%s %s]', sz_str, class(val));
end