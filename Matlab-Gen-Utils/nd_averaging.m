function [bin_values, bin_centers, bin_edges] = nd_averaging(val, coords, varargin)
%ND_AVERAGING Collapse a multi-dimensions data set into an array of average values
%
%   Given a set of data values defined in some multidimensional space, this
%   function will bin values along all the coordinate dimensions and, by
%   default, calculate the average of each bin.
%
%   [ BIN_VALUES, BIN_CENTERS, BIN_EDGES ] = ND_AVERAGING( VAL, COORDS )
%   Give the data values as VAL and the coordinates within the cell array
%   COORDS. Each entry in COORDS must be an array the same shape as VAL.
%
%   [ ___ ] = ND_AVERAGING( ___, BINS ) BINS defines what bins are to be
%   used; it must be a cell array with the same number of elements as
%   COORDS. Each element may define the bins for the corresponding
%   coordinate in one of three ways:
%
%       * If given as a scalar >= 1, that defines the number of bins in
%       that dimension. The bins will be evenly spaced between the smallest
%       and largest values.
%
%       * If given as a vector V with length N, then it will create N-1
%       bins, where the edges of bin N are V(N:N+1). V must be
%       monotonically increasing.
%
%       * If given as an N-by-2 matrix M, then it will define N bins with
%       edges M(N,:) for the Nth bin.
%
%   Example:
%
%       val = 1:8;
%       x = [1 1 1 1 2 2 2 2];
%       y = [10 10 20 20 10 10 20 20];
%       bins = {0.5:1:2.5, 5:10:20};
%       [bin_values, bin_centers] = nd_averaging(val, {x,y}, bins);
%
%   returns
%
%       bin_values = [  1.5 3.5  
%                       5.5 7.5  ]
%
%       bin_centers = { [1,2], [10,20] }
%
%
%   Parameters:
%
%       'op' - a string determining which operation to apply to the points
%       in each bin. Default is 'nanmean', i.e. all points in the bin are
%       averaged, excluding NaNs. May be any string that names a function
%       that accepts a vector of numbers and returns a scalar number.

p = advInputParser;
p.addOptional('bins', []);
p.addParameter('op','nanmean');

p.parse(varargin{:});
pout = p.Results;

bin_edges = pout.bins;
if isempty(bin_edges)
    bin_edges = make_default_bins(coords);
end

bin_op = setup_operation(pout.op);

bin_edges = setup_bins(bin_edges, coords);
output_dims = cellfun(@(x) size(x,1), bin_edges, 'uniform', false);
bin_values = nan(output_dims{:});

% Now we iterate over each bin, find the points in that bin, carry out the
% appropriate operation, and store the result in the corresponding
% location.

for i_bin = 1:numel(bin_values)
    % It's simpler to loop over a linear index than arbitrarily many
    % subscripts, so first we loop over the bins then since each dim of the
    % bin array corresponds to a different coordinate, we convert the
    % linear index to subscript indices so we can get the bin we're in for
    % each coordinate. 
    %
    % The use of a cell array in the output is just a way to get an
    % arbitrary number of outputs.
    bin_sub_inds = cell(size(bin_edges));
    [bin_sub_inds{:}] = ind2sub(size(bin_values), i_bin);
    
    % Now for each coordinate, we compare its values to the bin and
    % steadily remove points outside each bin in each dimension.
    xx = true(size(val));
    for i_coord = 1:numel(coords)
        this_bin = bin_edges{i_coord}(bin_sub_inds{i_coord},:);
        xx = xx & coords{i_coord} >= this_bin(1) & coords{i_coord} < this_bin(2);
    end
    
    % Then its just a matter of carrying out the averaging or other
    % operation.
    bin_values(i_bin) = bin_op(val(xx));
end

bin_centers = cellfun(@(x) mean(x,2), bin_edges, 'uniform', false);

end

function bins = make_default_bins(coords)
% By default, use 10 bins in each dimension
bins = repmat({10}, size(coords));
end

function bins = setup_bins(bins_in, coords)
E = JLLErrors;

if ~isequal(size(bins_in), size(coords)) || ~iscell(bins_in)
    E.badinput('If given, "bins" must be a cell array the same size as COORDS')
end

bins = cell(size(coords));
% Allow each coordinate to use a different type of bin definition: number
% of bins, vector, or matrix.
for i_coord = 1:numel(coords)
    this_bin = bins_in{i_coord};
    if isscalar(this_bin)
        if this_bin < 1
            E.badinput('Error in bin definition for coordinate #%d: scalar value < 1 given', i_coord)
        end
        bins{i_coord} = convert_num_bins_to_matrix(this_bin, coords{i_coord});
        
    elseif isvector(this_bin)
        if any(diff(this_bin) < 0)
            E.badinput('Error in bin definition for coordinate #%d: a vector must be monotonically increasing', i_coord);
        end
        bins{i_coord} = convert_vec_bins_to_matrix(this_bin);
        
    elseif ismatrix(this_bin) && size(this_bin,2) == 2
        bins{i_coord} = this_bin;
        
    else
        E.badinput('Error in bin definition for coordinate #%d: must be a scalar, vector, or n-by-2 matrix', i_coord);
    end
end

end

function bins = convert_num_bins_to_matrix(n_bins, coord_vals)
tmp_bins = linspace(min(coord_vals(:)), max(coord_vals(:))+eps, n_bins+1);
bins = convert_vec_bins_to_matrix(tmp_bins);
end

function bins = convert_vec_bins_to_matrix(bins_vec)
bins = [reshape(bins_vec(1:end-1),[],1), reshape(bins_vec(2:end),[],1)];
end

function op = setup_operation(op_name)
% If needed, we could check for certain custom strings that don't
% correspond to a Matlab function
op = str2func(op_name);

% str2func will return a handle to a function even if it doesn't exist.
% This will also catch if a function errors when given numbers. 
try
    chk = op(1:10);
catch err
    E = JLLErrors;
    E.badinput('The operation name "%s" gave the following error when the corresponding function was called with a vector of numbers: "%s"', op_name, err.message);
end

if ~isscalar(chk) || ~isnumeric(chk)
    E.badinput('The operation name "%s" did not return a scalar number', op_name);
end

end