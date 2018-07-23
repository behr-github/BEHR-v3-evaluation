function [ v ] = veccat( varargin )
%VECCAT Concatenate vectors along their non-singleton dimension
%   V = VECCAT( V1, V2, ... ) Concatenate the vectors V1, V2, etc. along
%   their long dimension. All inputs must be scalars or the same vector
%   type (i.e. row or column).
%
%   V = VECCAT( ___, 'column' ) forces the input vectors to be reshaped
%   into column vectors and concatenated along the first column.

E = JLLErrors;

xx = strcmpi(varargin, 'column');
if any(xx)
    force_column = true;
    varargin(xx) = [];
else
    force_column = false;
end

if ~force_column
    first_is_row = -1;
    for a=1:numel(varargin)
        % Checking if inputs are either all row or all column vectors. Cannot
        % just compare to first input, in case first input is scalar.
        if first_is_row < 0
            if isscalar(varargin{a}) || isempty(varargin{a})
                continue
            else
                first_is_row = isrow(varargin{a});
            end
        end
        if isrow(varargin{a}) ~= first_is_row && ~isscalar(varargin{a}) && ~isempty(varargin{a})
            E.badinput('All inputs to VECCAT must be the same sort of vector (row or column). This function does not handle mixed vector types')
        elseif ~isvector(varargin{a}) && ~isempty(varargin{a})
            E.badinput('VECCAT only concatenates vectors')
        end
    end
    
    non_empty_inputs = cellfun(@(x) ~isempty(x), varargin);
    
    % If all scalars, cat along first dimension
    if first_is_row < 0
        catdim = 1;
    else
        catdim = first_is_row + 1;
    end
    % Sometimes empty inputs behave as though their dimensions are not
    % consistent. Since concatenating them does nothing, just skip trying
    % to concatenate them.
    v = cat(catdim, varargin{non_empty_inputs});
else
    varargin = cellfun(@(x) reshape(x, [], 1), varargin, 'UniformOutput', false);
    v = cat(1, varargin{:});
end



end

