function [ gc_trop_mask ] = gc_trop_mask_by_pressure( gc_pres, gc_trop_pres )
%GC_TROP_MASK_BY_PRESSURE Returns a logical array where GC pressure >= tropopause pressure
%   Somewhat similar to set_stratosphere_nans, this function takes either
%   structures output from read_geos_output or their dataBlocks directly.
%   The first must contain the pressure levels of the GEOS-Chem boxes, the
%   second gives the tropopause pressure for each column.  It is required
%   that the first, second, and fourth dimensions of the first input
%   dataBlock be the same length as the first, second, and third dimensions
%   of the tropopause pressure dataBlock.  This will return a logical array
%   the same size as the pressure dataBlock.
%
%   Josh Laughner <joshlaugh5@gmail.com> 16 Nov 2015

E = JLLErrors;

if isstruct(gc_pres)
    if ~strcmp(gc_pres.fullName,'PEDGE-$_PSURF')
        E.badinput('First input must be a pressure level structure or the data matrix itself')
    end
    gc_pres = gc_pres.dataBlock;
end

if isstruct(gc_trop_pres)
    if ~strcmp(gc_trop_pres.fullName,'Tropopause pressure')
        E.badinput('First input must be a pressure level structure or the data matrix itself')
    end
    gc_trop_pres = gc_trop_pres.dataBlock;
end

if size(gc_pres,1) ~= size(gc_trop_pres,1) || size(gc_pres,2) ~= size(gc_trop_pres,2) || size(gc_pres,4) ~= size(gc_trop_pres,3)
    E.badinput('Lon, lat, and time dimensions of both input matrices must be the same')
end

gc_trop_mask = false(size(gc_pres));

for a=1:size(gc_pres,1)
    for b=1:size(gc_pres,2)
        for c=1:size(gc_pres,4)
            xx = squeeze(gc_pres(a,b,:,c)) >= gc_trop_pres(a,b,c);
            gc_trop_mask(a,b,xx,c) = true;
        end
    end
end
