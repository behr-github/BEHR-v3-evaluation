function [ NewMerge ] = combine_merges( varargin )
%COMBINE_MERGES Combines multiple "Merge" structures into one
%   NewMerge = COMBINE_MERGES( Merge1, Merge2, ... , MergeN ) Combines all
%   input Merge structures (with "metadata" and "Data" substructs) into
%   one, following these rules:
%       -- All Merges must represent the same date (as specified in the
%       metadata)
%       -- All Merges must have the same upper LoD and lower LoD flag
%       -- All UTC times that exist in the Merges will be represented
%       exactly once
%       -- Any field represented across multiple merges must have the same
%       fill value and unit

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

Merges = {};
param_args = {};

for a=1:numel(varargin)
    if isstruct(varargin{a}) && all(isfield(varargin{a}, {'Data', 'metadata'}))
        Merges{end+1} = varargin{a};
    else
        param_args{end+1} = varargin{a};
    end
end

p = inputParser;
p.addParameter('utc_field', 'UTC');
p.addParameter('new_utc_field', '');

p.parse(param_args{:});
pout = p.Results;
utc_field = pout.utc_field;
new_utc_field = pout.new_utc_field;

E = JLLErrors;

if ~ischar(utc_field)
    E.badinput('The value given for the parameter "utc_field" must be a string')
end

for a=1:numel(Merges)
    if ~isfield(Merges{a}.Data, utc_field)
        E.badinput('UTC field "%s" is not a Data field in merge #%d. Specify a field with the "utc_field" parameter if necessary.', a, utc_field)
    end
end

if ~ischar(new_utc_field)
    E.badinput('The value given for the parameter "new_utc_field" must be a string')
elseif isempty(new_utc_field)
    new_utc_field = utc_field;
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
E.addCustomError('bad_merge', 'combine_merge_fail', 'Failed combining merge: %s');

% First find all UTCs over all the Merges
utcs = [];
utc_ufn = get_unit_field(Merges{1}.Data.(utc_field));
for a=1:numel(Merges)
    if ~strcmp(Merges{a}.Data.(utc_field).(utc_ufn), Merges{1}.Data.(utc_field).(utc_ufn))
        E.callCustomError('bad_merge', 'Merge UTC fields are in different units');
    end
    utcs = veccat(utcs, Merges{a}.Data.(utc_field).Values);
end
utcs = sort(unique(utcs));

md_fns = fieldnames(Merges{1}.metadata);
md_struct = make_empty_struct_from_cell(md_fns, {{}});
NewMerge = struct('metadata', md_struct, 'Data', struct);

% Handle metadata. In most cases, just concatenate it from each merge,
% except for fields which really should be the same. An additional message
% explaining why the must be the same can be provided

req_same_val = {'upper_lod_flag', 'lower_lod_flag', 'date'};
additional_msg = {'', '', 'Merges spanning multiple dates are not supported.'};
for a=1:numel(Merges)
    for f=1:numel(md_fns)
        ff = strcmp(md_fns{f}, req_same_val);
        if any(ff)
            if a == 1
                NewMerge.metadata.(md_fns{f}) = Merges{a}.metadata.(md_fns{f});
            else
                if ~isequal(NewMerge.metadata.(md_fns{f}), Merges{a}.metadata.(md_fns{f}))
                    msg = sprintf('The metadata field "%s" has different values in different merges; this is not allowed for this field. %s', md_fns{f}, additional_msg{ff});
                    E.callCustomError('bad_merge', msg)
                end
            end
        else
            NewMerge.metadata.(md_fns{f}){end+1} = Merges{a}.metadata.(md_fns{f});
        end
    end
end

% Now handle data. Match by UTC value. Fix the incorrect capitalization of
% the UTC unit field
NewMerge.Data.(new_utc_field).Unit = Merges{1}.Data.(utc_field).(utc_ufn);
NewMerge.Data.(new_utc_field).Values = utcs;

utc_fills = cell(size(Merges));
for a=1:numel(Merges)
    utc_fills{a} = Merges{a}.Data.UTC.Fill;
end
utc_fills = unique(utc_fills);
if numel(utc_fills) > 1
    E.callCustomError('bad_merge', 'UTC fields have different fill values');
end
NewMerge.Data.(new_utc_field).Fill = utc_fills;

% Loop through each merge structure. Ignore UTC fields. Other fields,
% ensure that we are not overwriting existing values.

for a=1:numel(Merges)
    Data = Merges{a}.Data;
    fns = fieldnames(Data);
    for b=1:numel(fns)
        if strcmp(fns{b}, utc_field)
            continue
        elseif ~isfield(NewMerge.Data, fns{b});
            new_vals = Data.(fns{b}).Fill * ones(size(utcs));
            new_struct = struct('Unit', Data.(fns{b}).Unit, 'Fill', Data.(fns{b}).Fill, 'Values', new_vals);
            NewMerge.Data.(fns{b}) = new_struct;
        else
            % Verify that the unit and fill value are the same, if not we
            % cannot combine these structures b/c much of the existing code
            % relies on at least the Fill value being a scalar (so we can
            % just concatenate the different values into an array)
            
            if ~strcmp(NewMerge.Data.(fns{b}).Unit, Data.(fns{b}).Unit)
                msg = sprintf('Field %s has different units (%s) than in the combined merge (%s)', fns{b}, Data.(fns{b}).Unit, NewMerge.Data.(fns{b}).Unit);
                E.callCustomError('bad_merge', msg);
            elseif NewMerge.Data.(fns{b}).Fill ~= Data.(fns{b}).Fill
                msg = sprintf('Field %s has different fill value (%d) than in the combined merge (%d)', fns{b}, Data.(fns{b}).Fill, NewMerge.Data.(fns{b}).Fill);
                E.callCustomError('bad_merge', msg);
            end
        end
        
        % Now we just need to match values up by UTC
        for c=1:numel(Data.(fns{b}).Values)
            xx = find_new_utc(Data.(utc_field).Values(c), NewMerge.Data.(new_utc_field).Values, E);
            if NewMerge.Data.(fns{b}).Values(xx) ~= NewMerge.Data.(fns{b}).Fill
                msg = sprintf('Trying to assign a value for field %s at UTC %d, value already assigned', fns{b}, NewMerge.Data.(new_utc_field).Values(xx));
                E.callCustomError('bad_merge', msg);
            end
            NewMerge.Data.(fns{b}).Values(xx) = Data.(fns{b}).Values(c);
        end
    end
end
end

function ufield = get_unit_field(MField)
% Deal with broken capitalization for unit fields
mfns = fieldnames(MField);
rout = regexp(mfns, '[Uu]nit', 'once');
ufield = mfns{~iscellcontents(rout, 'isempty')};
end

function xx = find_new_utc(old_utc, new_utc, E)
if ~isscalar(old_utc)
    E.badinput('OLD_UTC must be a scalar')
end
xx = new_utc == old_utc;
if sum(xx) < 1
    msg = sprintf('Cannot find UTC %d from existing merge in combined merge.', old_utc);
    E.callCustomError('bad_merge', msg);
elseif sum(xx) > 1
    msg = sprintf('Multiple UTCs in the combined merge match UTC %d in the old merge', old_utc);
    E.callCustomError('bad_merge', msg);
end
end
