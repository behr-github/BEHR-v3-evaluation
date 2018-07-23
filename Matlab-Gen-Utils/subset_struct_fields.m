function S = subset_struct_fields(S, varargin)
%SUBSET_STRUCT_FIELDS Index all fields of a structure
%   S = SUBSET_STRUCT_FIELDS( S, LOGICAL1, LOGICAL2, ... ) iterates through
%   all fields in all elements of S and indexes that field with the first
%   LOGICAL that has the same size. Any number of LOGICAL arrays may be
%   given, but they must be logical arrays; linear or subscript indexing is
%   not supported.
%
%   Example:
%
%       s = struct('a',1:10,'b',rand(5,5),'c',rand(5,2),'d',reshape(1:25,5,5));
%       xx = 1:10 > 5;
%       s2 = subset_struct_fields(s, xx)
%       s2 = 
%
%         struct with fields:
%
%           a: [6 7 8 9 10]
%           b: [5x5 double]
%           c: [5x2 double]
%           d: [5x5 double]
%
%   Note that only the 'a' field was cut down because even though 'c' had
%   the same number of elements as xx, its shape was different.
%
%   Using the same s and xx:
%       yy = [false(5,3), true(5,2)];
%       s2 = subset_struct_fields(s, xx, yy)
%
%       s2 = 
%
%         struct with fields:
%
%           a: [6 7 8 9 10]
%           b: [10x1 double]
%           c: [5x2 double]
%           d: [10x1 double]
%
%   Both 'b' and 'd' were indexed by yy because it has the right shape.

fns = fieldnames(S);

for i_s = 1:numel(S)
    for i_field=1:numel(fns)
        for i_arg = 1:numel(varargin)
            logical_index = varargin{i_arg};
            sz_logical = size(logical_index);
            if isequal(size(S(i_s).(fns{i_field})), sz_logical)
                S(i_s).(fns{i_field}) = S(i_s).(fns{i_field})(logical_index);
                break
            end
        end
    end
end

end

