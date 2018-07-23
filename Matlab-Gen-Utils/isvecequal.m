function [ b ] = isvecequal( v1, v2 )
%ISVECEQUAL Tests if two vectors are equal as long as both are columns
%   B = ISVECEQUAL( V1, V2 ) Returns a scalar boolean, B, that is the
%   result of ISEQUAL(V1(:), V2(:)).


b = isequal(v1(:),v2(:));

end

