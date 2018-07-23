function [ frac ] = fracnan( M )
%FRACNAN Returns the fraction of values in M that are NaNs


frac = sum(isnan(M(:)))./numel(M);

end

