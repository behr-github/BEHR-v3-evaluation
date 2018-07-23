function sum_out = nanadd(varargin)
%NANADD Add arbritrarily many arrays, treating NaNs as zeros.
%   SUM = NANADD( V1, V2, V3, ... ) sum V1, V2, V3, etc. element-wise
%   treating NaNs as zeros unless an element is NaN in all inputs, in which
%   case it is kept a NaN. All inputs must be the same size.

E = JLLErrors;

target_size = size(varargin{1});
for i=2:numel(varargin)
    if ~isequal(size(varargin{i}), target_size)
        E.badinput('All inputs must have the same size');
    end
end

catdim = ndims(varargin{1}) + 1;

sum_out = nansum2(cat(catdim, varargin{:}), catdim);

end

