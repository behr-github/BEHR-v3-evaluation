function [ slon, slat ] = loadstates( subset )
%[ SLON, SLAT ] = LOADSTATES() Load in lines for US states.
%   [ SLON, SLAT ] = LOADSTATES returns lines for the continental states
%   alone by default.
%
%   [ SLON, SLAT ] = LOADSTATES('all') returns lines for all 50 states.
%
%   Josh Laughner <joshlaugh5@gmail.com> 18 Mar 2016

S = shaperead('usastatehi');

slon = [];
slat = [];

if ~exist('subset','var')
    subset = 'cont';
end

for a = 1:numel(S)
    if strcmpi(subset,'cont') && (strcmp(S(a).Name,'Alaska') || strcmp(S(a).Name,'Hawaii'))
        continue
    end
    slon = cat(2, slon, S(a).X);
    slat = cat(2, slat, S(a).Y);
end


end

