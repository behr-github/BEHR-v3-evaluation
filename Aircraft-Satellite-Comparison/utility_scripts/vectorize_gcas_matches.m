function data = vectorize_gcas_matches(data_in)
% VECTORIZE_GCAS_MATCHES Extract vectors of matched OMI and GCAS columns
%   DATA = VECTORIZE_GCAS_MATCHES( DATA_IN ) given a structure DATA_IN
%   returned by RUN_GCAS_VERIFICATION(), this will use the Matches field of
%   that structure to extract the pixels that contained GCAS data and
%   return a scalar structure DATA will all the fields of DATA_IN (except
%   matches) converted to vectors with all the matched data.

% We'll retain every field except the one that indicates whether a pixel
% had matching GCAS data.
E = JLLErrors;
fns = fieldnames(data_in);
fns(strcmp(fns, 'Matches')) = [];
data = make_empty_struct_from_cell(fns);

% Check all entries in data_in to determine which fields are 3D, if, for
% instance, a 20x1 chunk of pixels matches, then the pixel corner fields
% will end up being 2D
fields_dims = find_field_max_dims(data_in);
for a=1:numel(data_in)
    for b=1:numel(fns)
        val = data_in(a).(fns{b});
        if fields_dims.(fns{b}) == 2
            data.(fns{b}) = veccat(data.(fns{b}), val(data_in(a).Matches), 'column');
        elseif fields_dims.(fns{b}) == 3
            data.(fns{b}) = cat(2, data.(fns{b}), val(:, data_in(a).Matches));
        else
            E.notimplemented('Squeezing a %d dimensional variable not implemented', ndims(val));
        end   
    end
end
end

function fields_dims = find_field_max_dims(data)
fns = fieldnames(data);
fields_dims = make_empty_struct_from_cell(fns, 0);
for a=1:numel(data)
    for b=1:numel(fns)
        fields_dims.(fns{b}) = max(fields_dims.(fns{b}), ndims(data(a).(fns{b})));
    end
end
end