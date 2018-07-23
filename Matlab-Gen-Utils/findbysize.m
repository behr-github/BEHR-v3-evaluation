function [ indices ] = findbysize( X, n, varargin )
%findbysize Returns the indices of the smallest or largest n values in X.
%   Sorts X by value, keeping track of the initial linear indices of its
%   entries.  By default, finds the smallest n values in X, pass 'largest'
%   as the optional 3rd argument to reverse this behavior. NaNs are not
%   counted and are removed before finding the smallest or largest.


if numel(varargin)>0 && ischar(varargin{1});
    order = varargin{1};
elseif isempty(varargin);
    order = 'smallest';
else
    error(E.badinput('Third argument must be a string'));
end

Xvector = X(:); 
if(isrow(Xvector)); Xvector = Xvector'; end

M = [Xvector, (1:numel(Xvector))'];
M = sortrows(M);

nans = isnan(M(:,1));
M(nans,:) = [];

if strcmpi(order,'smallest');
    indices = M(1:n,2);
elseif strcmpi(order,'largest');
    indices = M(end-(n-1):end,2);
else
    error('findbysize:invalid_order','Third argument must be ''smallest'',''largest'', or empty');
end

% Return indices in ascending order as a row
indices = sortrows(indices)';

end

