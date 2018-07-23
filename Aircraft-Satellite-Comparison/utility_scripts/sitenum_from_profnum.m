function [ sitenums ] = sitenum_from_profnum( profnums )
%SITENUM_FROM_PROFNUM Convert DISCOVER-AQ profile numbers to site numbers
%   DISCOVER-AQ flights were constructed to give vertical profiles over
%   several sites in order to compare against satellite data. The profiles
%   were each assigned a unique ID number that contains the site number.
%   This number is given in the thousands place of that profile number;
%   this function extracts it.
%
%   SITENUMS = SITENUM_FROM_PROFNUM( PROFNUMS ) returns an array SITENUMS
%   the same size as PROFNUMS with the site numbers.

sitenums = floor(mod(profnums,100000)/1000);

end

