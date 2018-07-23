function S = cat_fields(varargin)
% CAT_FIELDS Concatenate vector fields in two or more structure
%   S = CAT_FIELDS(S1, S2, S3...) will concatenate the fields in all input
%   structures. All fields must be vectors, i.e. be able to be concatenated
%   with VECCAT. Example:
%
%       s1 = struct('a', 1:5, 'b', 10:10:50);
%       s2 = struct('a', 6:10, 'b', 60:10:100);
%
%       cat_fields(s1,s2)
%
%       ans = 
% 
%         struct with fields:
% 
%           a: [1 2 3 4 5 6 7 8 9 10]
%           b: [10 20 30 40 50 60 70 80 90 100]
%
%   By default, all structs must have the same fields and be the same size.
%
%   S = CAT_FIELDS( ___, 'strict' ) is the same as the default behavior
%   ('strict' requires that all fields are the same in the all structures).
%
%   S = CAT_FIELDS( ___, 'common' ) returns S which has only fields that
%   are common across all input structures. If no fields are common among
%   all input structures, an error is thrown.
%
%   S = CAT_FIELDS( ___, 'common-keep' ) is similar to 'common' except that
%   fields present in the first structure that are not common to all
%   structures are retained in S (but not concatenated in any way). This is
%   useful if you e.g. add some metadata fields to the first structure and
%   want those preseved.
%
%   S = CAT_FIELDS( ___, 'first' ) requires that all input structures have
%   all the fields present in the first structure. If not, an error is
%   thrown. Extra fields in later structures are allowed.

E = JLLErrors;

[all_structs, cat_fields, keep_fields] = parse_inputs(varargin);
S = all_structs{1};
S = remove_extra_fields(S, keep_fields);

for i_struct = 2:numel(all_structs)
    this_struct = all_structs{i_struct};
    for i_el = 1:numel(S)
        for i = 1:numel(cat_fields)
            try
                S(i_el).(cat_fields{i}) = veccat(S(i_el).(cat_fields{i}), this_struct(i_el).(cat_fields{i}));
            catch err
                if strcmp(err.identifier,'veccat:invalid_input')
                    % Give a clearer error message. This happens when we
                    % try to concatenate fields that aren't vectors.
                    E.notimplemented('%s', 'Concatenating non-vector fields not yet implemented');
                else
                    rethrow(err)
                end
            end
        end
    end
end
end

function [all_structs, cat_fields, keep_fields] = parse_inputs(inputs)
E = JLLErrors;

% Are all inputs (except maybe one) structs? Are at least two structs
% given?
xx = cellfun(@isstruct, inputs(:));
if sum(xx) < numel(inputs)-1 || (sum(xx) == numel(inputs)-1 && ~all(cellfun(@ischar, inputs(~xx))))
    E.badinput('At most one character input specifying the operation mode is allowed, all other inputs must be structs');
elseif sum(xx) < 2
    E.badinput('At least two structures required');
end

% Get the operating mode, if given. Otherwise default to 'strict';
allowed_op_modes = {'strict','common','common-keep','first'};
if sum(~xx) > 0
    op_mode = inputs{~xx};
    if ~ismember(op_mode, allowed_op_modes)
        E.badinput('The operating mode for CAT_FIELDS must be one of: %s', strjoin(allowed_op_modes, ', '));
    end
else
    op_mode = allowed_op_modes{1};
end

% Are all structs the same size?
all_structs = inputs(xx);
first_size = size(all_structs{1});
for i_struct = 2:numel(all_structs)
    this_size = size(all_structs{i_struct});
    if ~isequal(first_size, this_size)
        E.badinput('Size of structure #%d (%s) is not equal to that of the first structure (%s)', i_struct, this_size, first_size);
    end
end

% Now check the fields. If using 'strict' op mode, all structs must have
% the same fields. If using 'common', only copy the fields that are in
% common. If using 'first', then all fields in the first structure must be
% present in the later structs.
struct_fns = cellfun(@fieldnames, all_structs, 'UniformOutput', false);
switch op_mode
    case 'strict'
        [cat_fields, keep_fields] = check_fields_strict(struct_fns);
    case 'common'
        [cat_fields, keep_fields] = check_fields_common(struct_fns, false);
    case 'common-keep'
        [cat_fields, keep_fields] = check_fields_common(struct_fns, true);
    case 'first'
        [cat_fields, keep_fields] = check_fields_first(struct_fns);
    otherwise
        E.notimplemented('No action defined for op_mode = %s', op_mode);
end

end

function [cat_fns, keep_fns] = check_fields_strict(struct_fns)
% struct_fns should be a cell array of cell arrays, each sub-cell array is
% the field names of one of the input structs.
E = JLLErrors;

for i=2:numel(struct_fns)
    % Compare both directions (is fns_i in fns_1 and fns_1 in fns_i)
    % because we need to verify that all structures have the same fields.
    % Using ISMEMBER rather than STRCMP because we want to ignore field
    % order.
    missing_fields = any(~ismember(struct_fns{i}, struct_fns{1})) || any(~ismember(struct_fns{1}, struct_fns{i}));
    if missing_fields
        E.badinput('In "strict" operation mode, all given structs must have the same fields (though the order need not be the same)');
    end
end
cat_fns = struct_fns{1};
keep_fns = cat_fns;
end

function [cat_fns, keep_fns] = check_fields_common(struct_fns, keep_extras)
% struct_fns should be a cell array of cell arrays, each sub-cell array is
% the field names of one of the input structs.
E = JLLErrors;

xx = true(size(struct_fns{1}));
for i=2:numel(struct_fns)
    % Need to compare fns_1 in fns_i because we need to have ISMEMBER
    % output an array the same size as xx. We should only need to compare
    % one direction here, because we don't care if fns_i has extra fields.
    xx = xx & ismember(struct_fns{1}, struct_fns{i});
end

if sum(xx) == 0
    E.badinput('Input structs have no common fields');
end

cat_fns = struct_fns{1}(xx);
if keep_extras
    keep_fns = struct_fns{1};
else
    keep_fns = cat_fns;
end
end

function [cat_fns, keep_fns] = check_fields_first(struct_fns)
% struct_fns should be a cell array of cell arrays, each sub-cell array is
% the field names of one of the input structs.
E = JLLErrors;

for i=2:numel(struct_fns)
    % Comparing fns_1 is in fns_i because we need to verify that the fields
    % from the first structure are in all subsequent structures.
    xx = ~ismember(struct_fns{1}, struct_fns{i});
    if any(xx)
        E.badinput('In "first" operation mode, all given structs must at least contain the fields from the first struct');
    end
end

cat_fns = struct_fns{1};
keep_fns = cat_fns;
end

function S = remove_extra_fields(S, cat_fields)
s_fns = fieldnames(S);
xx = ~ismember(s_fns, cat_fields);
S = rmfield(S, s_fns(xx));
end
