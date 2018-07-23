function [ data, fill_value, scale_factor, offset ] = hdfreadmodis( filename, dsetname, varargin )
%HDFREADMODIS Wrapper around HDFREAD that handles fills, scale, and offset for MODIS files
%   DATA = HDFREADMODIS( FILENAME, DSETNAME ) - Reads in the dataset from
%   FILENAME at internal path DSETNAME. Uses HDFREAD internally as
%   HDFREAD(FILENAME, DSETNAME), so the DSETNAME must be formatted such
%   that HDFREAD understands it. Fill values defined by the _FillValue
%   attribute are converted to NaNs and the scale factor and offset are
%   applied as DATA = scale_factor * (hdf_values - add_offset).
%
%   DATA = HDFREADMODIS( ___, 'fill_crit', FILL_CRIT ) allows you to
%   specify the relative difference allowed between a value in the data set
%   and the fill value for the value to be considered a fill value, i.e.
%
%       abs((val - fill_val)/fill_val) < fill_crit 
%
%   must be met for VAL to be considered a fill value. Defaults to 0.001.
%
%   DATA = HDFREADMODIS( ___, 'log_index', {xx1, xx2, ...} ) uses the
%   logical vectors xx1, xx2, etc. to generate START, STRIDE, EDGE arrays
%   to subset the SDS during reading. xx1, xx2, etc. should be logical
%   arrays that are TRUE for indicies in the corresponding dimension of the
%   dataset that you want to include. Note: this assumes that the region to
%   read in is contiguous. If not, a warning is issued.
%
%   [ DATASET, FILL_VALUE, SCALE_FACTOR, OFFSET ] = HDFREADMODIS( __ )
%   returns the other attributes as well as the dataset.


E = JLLErrors;
p = inputParser;
p.addParameter('fill_crit',0.001);
p.addParameter('log_index', {});
p.parse(varargin{:});
pout = p.Results;

fill_crit = pout.fill_crit;
if ~isnumeric(fill_crit) || ~isscalar(fill_crit) || fill_crit <= 0
    E.badinput('FILL_CRIT must be a positive, scalar number')
end
log_index = pout.log_index;
if ~iscell(log_index) || any(~iscellcontents(log_index, 'isvector')) || any(~iscellcontents(log_index, 'islogical'))
    E.badinput('LOG_CRIT must be a cell array of logical vectors')
end

if isempty(log_index)
    data = double(hdfread(filename, dsetname));
else
    index_cell = make_index_cell(log_index);
    data = double(hdfread(filename, dsetname, 'index', index_cell));
end
fill_value = double(hdfreadatt(filename, dsetname, '_FillValue'));
try
    scale_factor = double(hdfreadatt(filename, dsetname, 'scale_factor'));
catch err
    if strcmp(err.identifier, 'hdfreadatt:hdf_attribute')
        warning('No scale_factor attribute found for %s, setting to 1', dsetname);
        scale_factor = 1;
    else
        rethrow(err)
    end
end

try
    offset = double(hdfreadatt(filename, dsetname, 'add_offset'));
catch err
    if strcmp(err.identifier, 'hdfreadatt:hdf_attribute')
        warning('No add_offset attribute found for %s, setting to 0', dsetname);
        offset = 0;
    else
        rethrow(err);
    end
end
fills = abs((data - fill_value) ./ fill_value) < fill_crit;
data(fills) = nan;

% This treatment is given by the MODIS MOD06 theoretical basis document, p.
% 87 (https://modis-atmos.gsfc.nasa.gov/_docs/C6MOD06OPUserGuide.pdf) and
% the metadata for MCD43C1 v006
% (https://ladsweb.modaps.eosdis.nasa.gov/api/v1/filespec/product=MCD43C1&collection=6,
% search BRDF_Albedo_Parameter).
data = (data - offset) * scale_factor;

end

function c = make_index_cell(log_index)
c = cell(1,3);
c([1,3]) = {zeros(1,numel(log_index))};
c(2) = {ones(1,numel(log_index))};

for a=1:numel(log_index)
    s = find(log_index{a}, 1, 'first');
    e = find(log_index{a}, 1, 'last')-s+1;
    if ~all(log_index{a}(s:e))
        warning('Subset of dimension %d is not contiguous, all elements between the first and last true value will be used', a)
    end
    c{1}(a) = s; % do not need to switch to 0 based indexing
    c{3}(a) = e;
end
end
