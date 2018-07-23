function [] = compare_profile(Merge,no2field,PROFILE,profnums,varargin)
%
%   Plots WRF-Chem NO2 profiles against in-situ aircraft profiles for
%   comparison.  Each will be normalized so that their maximum value = 1;
%   therefore we a really comparing shape factors.
%
%   Inputs:
%       Merge = an aircraft Merge data structure
%       no2field = a string indicating which field to use
%       PROFILE = a monthly WRF-Chem profile output with fields Latitude,
%           Longitude, Pressure, and NO2_Profile
%       profnums = profile numbers to use in the comparison.
%       [opt #1] = the name of the profile number field.
%
%   Josh Laughner <joshlaugh5@gmail.com> 5 Aug 2014

% Generate the grid matrices needed to interpolate the profile to the mean
% lat/lon of the aircraft profile

lonvec = PROFILE.Longitude(1,:);
latvec = PROFILE.Latitude(:,1);
pres = fliplr(PROFILE.Pressure);
no2 = flipdim(PROFILE.NO2_profile,1);
[Y,Z,X] = meshgrid(latvec,pres,lonvec);

% Load the data needed from the Merge structure

mlat = remove_merge_fills(Merge, 'LATITUDE');
mlon = remove_merge_fills(Merge, 'LONGITUDE');
mlon(mlon>180) = mlon(mlon>180)-360;
mno2 = remove_merge_fills(Merge, no2field);
mpres = remove_merge_fills(Merge, 'PRESSURE');

% Figure out if the user passed a profile field, if not, try to find one
if numel(varargin) > 0;
    proffield = varargin{1};
else
    proffields = search_merge_fields(Merge,'prof');
    if numel(proffields)==1; proffield = proffields{1};
    else error('compare_profile:profile_fields','Profile number field could not be identified.  Pass the field name as a string as the 5th argument');
    end
end

% Find all profile numbers that equal any of the profile numbers passed to
% this function
mprofnums = remove_merge_fills(Merge, proffield);
mm = false(size(mprofnums));
for a=1:numel(profnums)
    mm = mm | mprofnums == profnums(a);
end

% Restrict the data to the profiles in question
mlat = mlat(mm); mlon = mlon(mm);
mno2 = mno2(mm); mpres = mpres(mm);

% Interpolate the WRF profiles to the location of the aircraft profile. The
% order of the arguments may look funny, but it seems to work to match the
% ordering of the NO2 profile matrix.
compNO2prof = interp3(Y,Z,X,no2,nanmean(mlat),pres,nanmean(mlon));

% Bin the aircraft profiles to OMI pressure bins
[mno2_bins, mpres_bins] = bin_omisp_pressure(mpres, mno2);

% Normalize both NO2 profiles so the their max is 1
mno2_bins_norm = mno2_bins ./ max(mno2_bins);
compNO2prof_norm = compNO2prof ./ max(compNO2prof);

% Plot both profiles, the model as blue, the plane as red
figure;
l_model = line(compNO2prof_norm, pres, 'color', 'b', 'linewidth', 2, 'linestyle','--','marker','s');
l_airplane = line(mno2_bins_norm, mpres_bins, 'color', 'r', 'linewidth', 2, 'linestyle', ':', 'marker', '^');
set(gca,'ydir','reverse');
legend([l_model, l_airplane],{'WRF','In-situ'},'fontsize',16);
if numel(profnums) == 1;
    titleprofspec = '%d';
else
    titleprofspec = [repmat('%d, ',1,numel(profnums)-1), '%d'];
end

titlespec = ['Profiles ',titleprofspec,' on %s \nvs. WRF-Chem profiles'];
titlestr = sprintf(titlespec,profnums(:),Merge.metadata.date);
title(titlestr,'fontsize',16);
    
