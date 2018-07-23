function [ xx ] = containedin( A, B, case_sensitive )
%containedin Determines what elements of A are present in B.
%   Simple function that loops through each element in A and tests if 
%   any(A(i) == B).  Returns a logical matrix the same size as A. Can
%   accept matrices or cell arrays of strings. The option third argument is
%   a boolean for case sensitive comparison of cell arrays. Defaults to
%   true.

E = JLLErrors;

warning('The functionality of "containedin" is handled by the buildin MATLAB function "ismember". Please replace');

% Check input
if (~isnumeric(A) && ~iscell(A)) || (~isnumeric(B) && ~iscell(B))
    error(E.badinput('A and B can be matrices, vectors, scalars, or cell arrays only'));
end

if (isnumeric(A) && ~isnumeric(B)) || (iscell(A) && ~iscell(B))
    error(E.badinput('A and B must either both be matrices/vectors/scalars or both be cell arrays - they cannot be mixed'));
end

if iscell(A)
    for a=1:numel(A)
        if ~ischar(A{a})
            error(E.badinput('containedin only accepts cell arrays of strings'));
        end
    end
end

if iscell(B)
    for a=1:numel(B)
        if ~ischar(B{a})
            error(E.badinput('containedin only accepts cell arrays of strings'));
        end
    end
end

if nargin < 3
    case_sensitive = true;
elseif nargin >= 3
    if ~isscalar(case_sensitive)
        error(E.badinput('The optional third argument ''case_sensitive'' must be a scalar'));
    end
end

xx = false(size(A));
for i=1:numel(A)
    if isnumeric(A)
        if any(A(i) == B(:))
            xx(i) = true;
        end
    elseif iscell(A) && case_sensitive
        if any(strcmp(A{i},B))
            xx(i) = true;
        end
    elseif iscell(A) && ~case_sensitive
        if any(strcmpi(A{i},B))
            xx(i) = true;
        end
    else
        error(E.unknownError);
    end
end


end

