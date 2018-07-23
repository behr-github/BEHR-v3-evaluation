function [ dataBlock ] = set_stratosphere_nans( dataBlock, tp_level )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E=JLLErrors;

sz_data = size(dataBlock);
sz_data(3) = [];
sz_tp = size(tp_level);

if any(sz_data ~= sz_tp)
    E.badinput('lon, lat, and time dimensions of the data block and tp_level variables must match')
end

% Set a box that contains ANY stratospheric component to a NaN.
tp_level = floor(tp_level);

for t=1:size(tp_level,3);
    for a=1:size(tp_level,1);
        for b=1:size(tp_level,2);
            tp = tp_level(a,b,t);
            dataBlock(a,b,tp:end,t) = nan;
        end
    end
end

end

