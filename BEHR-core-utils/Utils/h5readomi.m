function [ dataset, fill_value, scale_factor, offset ] = h5readomi( filename, datasetname, varargin )
%H5READOMI Read a dataset in, removing fill values and applying scale and offset
%   DATASET = H5READOMI( FILENAME, DATASETNAME ) returns the dataset in
%   FILENAME specified by DATASETNAME with fill values converted into NaNs,
%   and the scale factor and offset applied as:
%       dataset = (dataset * scale_factor) - offset
%
%   DATASET = H5READOMI( FILENAME, DATASETNAME, 'fill_crit', FILL_CRIT )
%   allows you to specify the relative difference allowed between a value
%   in the data set and the fill value for the value to be considered a
%   fill value, i.e. 
%
%       abs((val - fill_val)/fill_val) < fill_crit 
%
%   must be met for VAL to be considered a fill value. Defaults to 0.001.
%
%   [ ___ ] = H5READOMI( ___, 'keep_type', true ) will keep the native type
%   of the value, rather than convert to a double. Note that integer types
%   will not have fill values replaced with NaNs, because they do not
%   support NaNs.
%
%   [ DATASET, FILL_VALUE, SCALE_FACTOR, OFFSET ] = H5READOMI( __ ) returns
%   the other attributes as well as the dataset.

E = JLLErrors;
p = inputParser;
p.addParameter('fill_crit',0.001);
p.addParameter('keep_type', false);
p.parse(varargin{:});
pout = p.Results;

fill_crit = pout.fill_crit;
if ~isnumeric(fill_crit) || ~isscalar(fill_crit) || fill_crit <= 0
    E.badinput('FILL_CRIT must be a positive, scalar number')
end

keep_type = pout.keep_type;
if ~islogical(keep_type) || ~isscalar(keep_type)
    E.badinput('KEEP_TYPE muts be a scalar logical')
end

dataset = h5read(filename, datasetname);
fill_value = h5readatt(filename, datasetname, '_FillValue');
scale_factor = h5readatt(filename, datasetname, 'ScaleFactor');
offset = h5readatt(filename, datasetname, 'Offset');

if ~keep_type
    dataset = double(dataset);
    fill_value = double(fill_value);
    scale_factor = double(scale_factor);
    offset = double(offset);
end

if ~isinteger(dataset)
    fills = abs((dataset - fill_value) ./ fill_value) < fill_crit;
    dataset(fills) = nan;
end

dataset = (dataset*scale_factor) - offset;

% I can't find any source that explicitly indicates whether scale or offset
% should be applied first. This order is how its done in the PSM code. The
% following warning will be issued if both matter.

if scale_factor ~= 1 && offset ~= 0
    warning('Dataset %s in file %s has both scaling and offset to apply, the order of operations is currently scale, then offset, but this is unconfirmed',...
        datasetname, filename)
end
end
