function avg = easymean(varargin)
%EASYMEAN Quickly average multiple arrays of the same size
%   AVG = EASYMEAN( M1, M2, ... ) Given multiple arrays M1, M2, etc. that
%   are the same size, this will average them using NANMEAN(). AVG will be
%   an array the same size as all the inputs.

E = JLLErrors;

if numel(varargin) < 1
    E.badinput('Must give at least one input')
end

all_sizes = cellfun(@size, varargin, 'uniform', false);
if ~all(cellfun(@(x) isequal(x, all_sizes{1}), all_sizes))
    E.badinput('All given arrays must have the same shape');
end

cat_dim = ndims(varargin{1})+1;
avg = nanmean(cat(cat_dim, varargin{:}), cat_dim);

end

