function [ dataBlock ] = set_lightning_nans( dataBlock, lnox_logical )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;

for t=1:size(lnox_logical,3)
    crit = lnox_logical(:,:,t);
    for l=1:size(dataBlock,3)
        level = dataBlock(:,:,l,t);
        level(~crit) = nan;
        dataBlock(:,:,l,t) = level;
    end
end
end

