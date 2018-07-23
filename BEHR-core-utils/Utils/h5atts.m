function [ attrs ] = h5atts( varargin )
%H5ATTS Lists dataset attributes from and HDF5 file
%   ATTRS = H5ATTS( INFO, GRP_IND, GRP_IND, ..., DATASET ) Returns the
%   attributes names as values as strings in the cell array ATTRS for a
%   dataset in the file represented by the INFO structure, which must be
%   one returned from H5INFO(). The following arguments, GRP_IND, are the
%   series of group indices in the hierarchy of the HDF5 file. DATASET is
%   the index of the dataset (or group) to list attributes for. For
%   example, if the dataset of interest can be reached with:
%
%       info.Groups(1).Groups(2).Groups(1).Datasets(10)
%
%   call this function as:
%
%       h5atts(info, 1, 2, 1, 10)
%
%   ATTRS = H5ATTS( ___, 'unannotated' ) returns the attribute information
%   instead as a struct, where each field name is an attribute name and the
%   value the attribute value.

E = JLLErrors;

xx = strcmpi('unannotated', varargin);
if sum(xx) > 0
    unannot = true;
    varargin(xx) = [];
else
    unannot = false;
end

if nargin == 1 %If there's only one argument, assume that it is the dataset we want to enumerate the attributes for
    h5dataset = varargin{1};
    
else
    h5group = varargin{1}; %Save the h5 object (info or group)
    n = length(varargin) - 2;
    groups = cell2mat(varargin(2:end));
    for i=1:n
        h5group = h5group.Groups(groups(i));
    end
    
    i = varargin{end};
    % Assume that the last index is pointing to a dataset by default. If
    % there aren't any datasets, then it's probably pointing to a group.
    if isempty(h5group.Datasets)
        h5dataset = h5group.Groups(i);
    else
        h5dataset = h5group.Datasets(varargin{end});
    end
end


n = length(h5dataset.Attributes);
if ~unannot
    attrs = cell(n+1,1);
    attrs{1} = sprintf('Dataset: %s', h5dataset.Name);
else
    attrs = make_empty_struct_from_cell({h5dataset.Attributes.Name});
    fns = fieldnames(attrs); % some field names will have been changed to work as structure fields
end
for i = 1:n
    if ~unannot
        if isnumeric(h5dataset.Attributes(i).Value)
            val = mat2str(h5dataset.Attributes(i).Value);
        elseif ischar(h5dataset.Attributes(i).Value)
            val = h5dataset.Attributes(i).Value;
        elseif iscellstr(h5dataset.Attributes(i).Value)
            val = strjoin(h5dataset.Attributes(i).Value, ', ');
        else
            E.notimplemented('Attribute value is not numeric, character, or cellstr');
        end
        info = sprintf('%d: %s = %s', i, h5dataset.Attributes(i).Name, val);
        attrs{i+1} = info;
    else
        attrs.(fns{i}) = h5dataset.Attributes(i).Value;
    end
    
end

end

