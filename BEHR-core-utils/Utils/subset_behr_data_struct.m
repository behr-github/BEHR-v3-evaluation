function Data = subset_behr_data_struct(Data, varargin)
% SUBSET_BEHR_DATA_STRUCT Subset fields of a BEHR data structure, keeping shape as much as possible
%   DATA = SUBSET_BEHR_DATA_STRUCT(DATA, INDEX1, INDEX2, ...) Subsets the
%   fields in DATA with INDEX1, INDEX2, etc. Each INDEX argument
%   corresponds to one element in DATA, so if DATA is a scalar struct, only
%   one INDEX argument is needed, if DATA is a 4x1 struct, four INDEX
%   arguments are needed. Each index must be either a logical or linear
%   index array for one of the 2D fields in DATA. 
%
%   Each field is subset according to the following rules:
%       * A 2D field is directly indexed with the INDEX argument for its
%       element of DATA.
%
%       * A 1D field with along- or across- track dimensions is expanded
%       into a 2D field before indexing with the 2D index array.
%
%       * A 3D field 

E = JLLErrors;
% First cut down the various 2-, and 3-D fields. There are both along- and
% across- track 1D fields, so handle them by replicating them into 2D
% matrices, and there are 3D fields where the first dimension is either 4
% (corners) or ~30 (scattering weights, etc).

if ~isstruct(Data)
    E.badinput('DATA must be a structure')
end

if numel(varargin) ~= numel(Data)
    E.badinput('Number of indexing arrays must match the number of elements in Data')
end

fns = fieldnames(Data);
for i_data = 1:numel(Data)
    this_data = Data(i_data);
    size_2d = size(this_data.Longitude);
    if size_2d(1) == size_2d(2)
        E.notimplemented('%s', 'Cannot handle an orbit with equal along- and across- track dimensions');
    end
    
    n_corners = size(this_data.FoV75CornerLongitude,1);
    n_plevs = size(this_data.BEHRPressureLevels,1);
    
    [xx_2d, xx_corners, xx_plevs] = format_input_indexing_array(size_2d, varargin{i_data}, n_corners, n_plevs, i_data);
    
    leading_dim_size = zeros(size(fns));
    
    for i_fn = 1:numel(fns)
        if isvector(this_data.(fns{i_fn}))
            % Assume that along-track only arrays are already column
            % vectors and across-track only arrays
            if numel(this_data.(fns{i_fn})) == size_2d(1)
                this_data.(fns{i_fn}) = repmat(this_data.(fns{i_fn}),1,size_2d(2));
            elseif numel(this_data.(fns{i_fn})) == size_2d(2)
                this_data.(fns{i_fn}) = repmat(this_data.(fns{i_fn}),size_2d(1),1);
            end
            % If the 1D array does not match the along or across track
            % dimension, do nothing.
        elseif ndims(this_data.(fns{i_fn})) == 3
            % 3D fields with a first dimension != corners or plevs will not
            % be subset, so when they get reshaped at the end, they'll just
            % end up being 2D with the last two dims rearranged into one.
            leading_dim_size(i_fn) = size(this_data.(fns{i_fn}),1);
        end
    end
    
    this_data = subset_struct_fields(this_data, xx_2d, xx_corners, xx_plevs);
    for i_fn = 1:numel(fns)
        % Reshape 3D fields so that their first dimension is the same as it
        % was originally, but the next two dimensions are compressed to
        % one.
        if leading_dim_size(i_fn) > 0
            this_data.(fns{i_fn}) = reshape(this_data.(fns{i_fn}), leading_dim_size(i_fn), []);
        end
    end
    
    Data(i_data) = this_data;
end

end

function [xx, xx_corners, xx_plevs] = format_input_indexing_array(size_2d, input_index, n_corners, n_plevs, i_data)
% i_data just used to make the error messages more intelligent
E = JLLErrors;
if islogical(input_index)
    if ~isequal(size(input_index), size_2d)
        E.badinput('The %1$s input index is a different size than the 2D arrays in the %1$s element of Data', ordinal_string(i_data));
    else
        xx = input_index;
    end
elseif isnumeric(input_index)
    xx = false(size_2d);
    xx(input_index) = true;
else
    E.badinput('The %s input index is neither a logical nor numeric array', ordinal_string(i_data));
end

xx_corners = repmat(permute(xx,[3 1 2]),n_corners,1,1);
xx_plevs = repmat(permute(xx,[3 1 2]),n_plevs,1,1);

end

function s = ordinal_string(i)
% Return an ordinal number string. Make it so that e.g. 1 = 1st, 101 =
% 101st, 2 = 2nd, 102 = 102nd, etc.
if mod(abs(i),100) == 1
    suffix = 'st';
elseif mod(abs(i),100) == 2
    suffix = 'nd';
elseif mod(abs(i),100) == 3
    suffix = 'rd';
else
    suffix = 'th';
end

s = sprintf('%d%s',i,suffix);

end