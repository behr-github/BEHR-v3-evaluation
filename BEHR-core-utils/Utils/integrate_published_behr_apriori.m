function [ vcds ] = integrate_published_behr_apriori( no2_profs, pres_levs, surf_pres )
%INTEGRATE_PUBLISHED_BEHR_APRIORI Calculate VCDs from published a priori
%   VCDS = INTEGRATE_PUBLISHED_BEHR_APRIORI( DATA ) calculates VCDs for the
%   BEHRNO2apriori field in the structure DATA.
%
%   VCDS = INTEGRATE_PUBLISHED_BEHR_APRIORI( NO2_PROFS, PRES_LEVS,
%   SURF_PRES ) calculates VCDs for NO2 profiles NO2_PROFS (in mixing
%   ratio) defined as pressures PRES_LEVS with surface pressures SURF_PRES
%   (both pressure fields in hPa). NO2_PROFS and PRES_LEVS must be 3D
%   arrays with the vertical dimension in the first array and the same
%   overall shape. SURF_PRES must be a 2D array the same size as the second
%   and third dimensions of NO2_PROFS and PRES_LEVS.
%
%   This function is used to integrate the NO2 a priori profiles we publish
%   as part of the BEHR product. It uses integPr2 internally, so the
%   integration is the same as in the main BEHR algorithm. The only tricky
%   part is handling NaNs in the pressure and NO2 vectors; any indices with
%   a NaN in either are removed and if there are not at least 2 non-NaN
%   values, the integration will not be carried out (a NaN is given in the
%   VCDs).

E = JLLErrors;

if nargin == 1
    Data = no2_profs;
    req_fields = {'BEHRNO2apriori', 'BEHRPressureLevels', 'GLOBETerpres'};
    if ~isstruct(Data) || ~isscalar(Data)
        E.badinput('DATA must be a scalar structure (in the one argument form)');
    elseif any(~isfield(Data, req_fields))
        E.badinput('DATA must have the fields: %s', strjoin(req_fields, ', '));
    end
    
    no2_profs = Data.BEHRNO2apriori;
    pres_levs = Data.BEHRPressureLevels;
    surf_pres = Data.GLOBETerpres;
else
    if ~isnumeric(no2_profs)
        E.badinput('NO2_PROFS must be a numeric array');
    elseif ~isnumeric(pres_levs)
        E.badinput('PRES_LEVS must be a numeric array');
    elseif ~isnumeric(surf_pres)
        E.badinput('SURF_PRES must be a numeric array');
    end
    
    sz_profs = size(no2_profs);
    sz_preslevs = size(pres_levs);
    sz_surfpres = size(surf_pres);
    
    if ~isequal(sz_profs, sz_preslevs)
        E.badinput('NO2_PROFS and PRES_LEVS must be the same size')
    elseif ~isequal(sz_profs(2:end), sz_surfpres)
        E.badinput('SURF_PRES must be the same size as the last two dimensions of NO2_PROFS')
    end
end

vcds = nan(size(surf_pres));
for a=1:numel(surf_pres)
    notnans = ~isnan(no2_profs(:,a)) & ~isnan(pres_levs(:,a));
    if sum(notnans) < 2
        continue
    end
    
    vcds(a) = integPr2(no2_profs(notnans,a), pres_levs(notnans,a), surf_pres(a));
end

end

