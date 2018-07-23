function [ xxa, xxb ] = find_common_elements( A, B, varargin )
%FIND_COMMON_ELEMENTS Returns indices of common elements in A and B
%   [ XXA, XXB ] = FIND_COMMON_ELEMENTS( A, B ) Returns vectors XXA and XXB
%   which are indicies for A and B respectively that are elements found in
%   both A and B, that is, A(XXA) will return elements that exist in B and
%   B(XXB) will return elements that exist in A. A and B can be numeric or
%   cell arrays. The indicies of XXA will always be in order, the indicies
%   of XXB will be in the order that those elements appear in A.
%
%   [ XXA, XXB ] = FIND_COMMON_ELEMENTS( ___, 'nodup' ) does not include
%   duplicated elements in the indices, instead only returning the first
%   occurance of the common element. That is, for A = [10 20 20 40] and B =
%   [10 20], XXA will only be [1, 2]. This also ensures that A(XXA) and
%   B(XXB) will be the same array.
%
%   [ XXA, XXB ] = FIND_COMMON_ELEMENTS( ___, 'ignorecase' ) causes any
%   string comparisons in cell arrays to be case insensitive.

nodup=false;
nocase=false;
if ismember('nodup',varargin)
    nodup = true;
end
if ismember('ignorecase',varargin)
    nocase = true;
end

E = JLLErrors;
if (~iscell(A) && ~isnumeric(A)) || (~iscell(B) && ~isnumeric(B))
    E.badinput('A and B should be numeric or cell arrays')
elseif xor(iscell(A), iscell(B))
    E.badinput('A and B must both be cell arrays OR both be numeric arrays, not mixed')
end


if isnumeric(A) && isnumeric(B)
    comp_fxn = @(A,B,i) compare_num(A,B,i);
elseif iscellstr(A) && iscellstr(B)
    comp_fxn = @(A,B,i) compare_cellstr(A,B,i);
else
    comp_fxn = @(A,B,i) compare_cell(A,B,i);
end

xxa = [];
xxb = [];
for i=1:numel(A)
    xx_i = comp_fxn(A,B,i);
    if sum(xx_i) > 0
        xxa = cat(2, xxa, i);
        if nodup
            xxb = cat(2, xxb, find(xx_i,1,'first'));
        else
            xxb = cat(2, xxb, find(xx_i));
        end
    end
end

    
    function xx = compare_num(A, B, ind)
        xx = A(ind) == B;
    end

    function xx = compare_cellstr(A, B, ind)
        if nocase
            xx = strcmpi(A{ind}, B);
        else
            xx = strcmp(A{ind}, B);
        end
    end

    function xx = compare_cell(A,B,ind)
        xx = false(size(B));
        if nocase && ischar(A{ind})
            a_i = lower(A{ind});
        else
            a_i = A{ind};
        end
        for j=1:numel(B)
            if nocase && ischar(B{j})
                b_j = lower(B{j});
            else
                b_j = B{j};
            end
            if isequal(a_i,b_j)
                xx(j) = true;
            end
        end
    end

end

