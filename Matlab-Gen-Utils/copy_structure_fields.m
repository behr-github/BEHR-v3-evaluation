function [ dest ] = copy_structure_fields(source, dest, varargin)
%COPY_STRUCTURE_FIELDS Copies fields values from one struct to another
%   DEST_OUT = COPY_STRUCTURE_FIELDS( SOURCE, DEST ) Given two structures,
%   SOURCE and DEST, return a structure with the same fields as DEST where
%   each field has the value of that field in SOURCE.
%
%   DEST_OUT = COPY_STRUCTURE_FIELDS( SOURCE, FIELDS_TO_COPY ) If instead
%   given a cell array of strings FIELDS_TO_COPY, those fields will be
%   copied from SOURCE into a structure DEST_OUT with those fields.
%
%   DEST_OUT = COPY_STRUCTURE_FIELDS( SOURCE, DEST, FIELDS_TO_COPY ) This
%   form returns DEST_OUT with all the same fields as the input DEST
%   structure but only copies the fields specified in the cell array
%   FIELDS_TO_COPY. If there are additional fields in DEST not included in
%   FIELDS_TO_COPY, they retain their input values in DEST_OUT. (This is
%   useful if trying to place the output in a non-scalar structure and so
%   it must have the same fields as the rest of the structure, but not all
%   those fields are in SOURCE or you do not want to overwrite all of the
%   fields.)
%
%   DEST_OUT = COPY_STRUCTURE_FIELDS( ___, 'substructs' ) will search not
%   just SOURCE for the fields to copy but also any substructures within
%   SOURCE. I.e., the fields to copy do not need to be in the top level of
%   SOURCE. This may be used with any of the above syntaxes.
%
%   DEST_OUT = COPY_STRUCTURE_FIELDS( ___, 'ignore_error' ) will ignore
%   errors caused by trying to copy a field from SOURCE that does not
%   exist. Instead, the original value (if DEST is given as a structure in
%   the input) or the default value of an empty array (if just
%   FIELDS_TO_COPY is given) is retained for that field.
%
%   DEST_OUT = COPY_STRUCTURE_FIELDS( ___, 'missing' ) will copy fields
%   from SOURCE that are missing from DEST. This ignores FIELDS_TO_COPY in
%   the three argument form and errors if used with the two argument form.

%%%%%%%%%%%%%%%%%
% INPUT PARSING %
%%%%%%%%%%%%%%%%%

E = JLLErrors;

p = advInputParser;
p.addOptional('fields_to_copy', {}, @iscellstr);
p.addFlag('substructs');
p.addFlag('ignore_error');
p.addFlag('missing');

p.parse(varargin{:});
pout = p.AdvResults;

ignore_error = pout.ignore_error;
missing_flag = pout.missing;

if ~isstruct(source)
    E.badinput('SOURCE must be a structure');
end

if iscellstr(dest) || ischar(dest)
    if missing_flag
        E.badinput('With the "missing" flag, DEST must be a struct');
    elseif ischar(dest)
        dest = {dest};
    end
    fields_to_copy = dest;
    dest = make_empty_struct_from_cell(dest, fields_to_copy);
    dest = repmat(dest, size(source));
elseif isstruct(dest)
    if ~isequal(size(dest), size(source))
        E.badinput('DEST must the same size as SOURCE, if given as a structure');
    end
    
    if ~missing_flag
        fields_to_copy = pout.fields_to_copy;
        % If fields_to_copy isn't specified then we should copy all fields
        % present in the destination structure
        if isempty(fields_to_copy)
            fields_to_copy = fieldnames(dest);
        end
    else
        fields_to_copy = missing_fields(source, dest);
        % If only copying missing fields, then if none are missing, there's
        % nothing to copy.
    end
    
else
    E.badinput('DEST must be a string, cell array of strings, or structure');
end

if pout.substructs
    finder_fxn = @find_substruct_field;
    err_addendum = ' or any of its substructures';
else
    finder_fxn = @find_field;
    err_addendum = '';
end

%%%%%%%%%%%%%%%%%
% MAIN FUNCTION %
%%%%%%%%%%%%%%%%%

for i_ind = 1:numel(dest)
    for i_field = 1:numel(fields_to_copy)
        try
            dest(i_ind).(fields_to_copy{i_field}) = finder_fxn(source(i_ind), fields_to_copy{i_field});
        catch err
            if strcmp(err.identifier, 'MATLAB:nonExistentField')
                % Reformat the message to better reflect the true cause of the
                % error
                if ~ignore_error
                    error('MATLAB:nonExistentField', 'Could not find field "%s" in SOURCE%s', fields_to_copy{i_field}, err_addendum);
                end
            else
                rethrow(err)
            end
        end
    end
end


end

function fns = missing_fields(source, dest)
% Returns a cell array of fields in SOURCE but not DEST
fns_source = fieldnames(source);
fns_dest = fieldnames(dest);
xx = ~ismember(fns_source, fns_dest);
fns = fns_source(xx);
end

function val = find_field(S, field_name)
val = S.(field_name);
end