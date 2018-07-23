function [ Matched_Data ] = gcas_sao_verification( gcas_files, Data_in, varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

E = JLLErrors;
p = inputParser;
p.addParameter('cloud_prod', 'omi');
p.addParameter('cloud_frac_max', 0.2);
p.addParameter('row_anomaly', 'XTrackFlags');
p.addParameter('sat_fields', {'BEHRColumnAmountNO2Trop', 'ColumnAmountNO2Trop'});
p.addParameter('time_window', 1.5);
p.addParameter('DEBUG_LEVEL', 2);

p.parse(varargin{:});
pout = p.Results;

cloud_prod = pout.cloud_prod;
cloud_frac_max = pout.cloud_frac_max;
row_anomaly = pout.row_anomaly;
sat_fields = pout.sat_fields;
time_window = pout.time_window;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

gcas_files = files_input(gcas_files);
GCAS = read_gcas_files(gcas_files);

for d=1:numel(Data_in)
    % Reject satellite
    Data_in(d).Areaweight = ones(size(Data_in(d).Longitude));
    if isfield(Data_in, 'BEHRQualityFlags')
        % The BEHRQualityFlags field was added in v3.0A of BEHR, so as long
        % as it's present, we can use the newest version of
        % omi_pixel_reject.
        if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject.m\n'); end
        reject_detail = struct('cloud_type', cloud_prod, 'cloud_frac', cloud_frac_max,...
            'row_anom_mode', row_anomaly, 'check_behr_amf', true);
        this_data = omi_pixel_reject(Data_in(d),'detailed',reject_detail);
    elseif isfield(Data_in, 'BEHRColumnAmountNO2Trop')
        % If there's a BEHR NO2 column field but not the BEHR quality flags
        % field, we're probably working with version 2 data.
        if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_pixel_reject_v2.m\n'); end
        this_data = omi_pixel_reject_v2(Data_in(d), cloud_prod, cloud_frac_max, row_anomaly);
    else
        % If the BEHR NO2 column field isn't in Data, it's definitely not a
        % BEHR file, so use the SP reject field.
        if DEBUG_LEVEL > 0; fprintf('   Rejecting with omi_sp_pixel_reject.m\n'); end
        this_data = omi_sp_pixel_reject(Data_in(d),cloud_prod,cloud_frac_max,row_anomaly);
    end
    
    xx_reject = this_data.Areaweight == 0;
    for a=1:numel(sat_fields)
        this_data.(sat_fields{a})(xx_reject) = nan;
    end
    
    % Cut down the satellite orbit to just the area that is relevant for
    % the GCAS data
    this_data = cut_down_satellite_orbit(this_data, GCAS, sat_fields);
    
    % Match the GCAS data in space and time
    Matched_Data(d) = match_gcas_to_omi(this_data, GCAS, time_window);
end

end

function [GCAS] = read_gcas_files(gcas_files)
GCAS.Longitude = [];
GCAS.Latitude = [];
GCAS.NO2vcd = [];
GCAS.CloudFlag = [];
GCAS.UTCDateTime = [];

% This structure should have the same fields as GCAS each of which points
% to a cell array, C, that if given to h5dsetname as h5dsetname(hi, C{:})
% will return the right dataset. Fields in GCAS but not here must be
% handled manually
var_mapping = struct('Longitude', {{2, 'Longitude'}},...
                     'Latitude', {{2, 'Latitude'}},...
                     'NO2vcd', {{1, 'NO2_RetrievedVerticalColumnBelow'}},...
                     'CloudFlag', {{1, 'CloudFlag'}});

fns = fieldnames(var_mapping);

                 
for a=1:numel(gcas_files)
    hi = h5info(gcas_files{a});
    % There does not appear to be any scaling or offset in these data
    for f=1:numel(fns)
        dset_name = h5dsetname(hi, var_mapping.(fns{f}){:});
        %dset = h5info(hi.Filename, dset_name);
        val = double(h5read(hi.Filename, dset_name));
        % Given that the fill values are 0, they don't seem to actually use
        % them
        %val(val == double(dset.FillValue)) = NaN;
        GCAS.(fns{f}) = cat(1, GCAS.(fns{f}), val(:));  
    end
    
    % Convert the cloud flag to a logical
    GCAS.CloudFlag = logical(GCAS.CloudFlag);
    
    % Convert GCAS date and UTC time into a Matlab datenum
    yr = double(h5read(hi.Filename, h5dsetname(hi,2,'Year')));
    mn = double(h5read(hi.Filename, h5dsetname(hi,2,'Month')));
    dy = double(h5read(hi.Filename, h5dsetname(hi,2,'Day')));
    utc_hr_of_day = double(h5read(hi.Filename, h5dsetname(hi,2,'UTCHours')));
    
    utc_datetime = datenum(yr,mn,dy) + utc_hr_of_day/24;
    GCAS.UTCDateTime = cat(1, GCAS.UTCDateTime, utc_datetime(:));
    
end

% Go ahead and set the NO2 to NaN where the retrieval is over clouds
GCAS.NO2vcd(GCAS.CloudFlag) = nan;

end

function data_out = cut_down_satellite_orbit(data, GCAS, sat_fields)
sat_loncorn = reshape(data.FoV75CornerLongitude,4,[]);
sat_latcorn = reshape(data.FoV75CornerLatitude,4,[]);

gcas_lon_bounds = [min(GCAS.Longitude(:)), max(GCAS.Longitude(:))];
gcas_lat_bounds = [min(GCAS.Latitude(:)), max(GCAS.Latitude(:))];

gcas_xv = [gcas_lon_bounds(1), gcas_lon_bounds(1), gcas_lon_bounds(2), gcas_lon_bounds(2)];
gcas_yv = [gcas_lat_bounds(1), gcas_lat_bounds(2), gcas_lat_bounds(2), gcas_lat_bounds(1)];

xx_in = any(inpolygon(sat_loncorn, sat_latcorn, gcas_xv, gcas_yv),1);
xx_in = reshape(xx_in, size(data.Longitude));
xx = any(xx_in, 2);
yy = any(xx_in, 1);

% Specify additional 2D fields that should be copied here
copy_fields = cat(1, {'Longitude';'Latitude'}, sat_fields(:));

% 1- or 3- D fields that need special handling. Go ahead and convert OMI
% time from TAI93 to a Matlab datenum.
time = repmat(data.Time, 1, size(data.Longitude,2));
data_out.UTCDateTime = omi_time_conv(time(xx,yy));
data_out.FoV75CornerLongitude = data.FoV75CornerLongitude(:,xx,yy);
data_out.FoV75CornerLatitude = data.FoV75CornerLatitude(:,xx,yy);
for a=1:numel(copy_fields)
    data_out.(copy_fields{a}) = data.(copy_fields{a})(xx,yy);
end


end

function [this_data] = match_gcas_to_omi(this_data, GCAS, time_window_hours)
% First write the GCAS fields into this_data, to ensure that we always
% return a structure with those field
fns = fieldnames(GCAS);
default_mat = nan(size(this_data.Longitude));
for f=1:numel(fns)
    new_fn = sprintf('GCAS_%s',fns{f});
    this_data.(new_fn) = default_mat;
end
this_data.Matches = false(size(this_data.Longitude));

if isempty(this_data.Longitude)
    return
end
% Next cut down the GCAS data to the time range that we specified in the
% input
tt = GCAS.UTCDateTime >= (min(this_data.UTCDateTime(:)) - time_window_hours/24) & GCAS.UTCDateTime <= (max(this_data.UTCDateTime(:)) + time_window_hours/24);

for f=1:numel(fns)
    GCAS.(fns{f}) = GCAS.(fns{f})(tt);    
end

% Now loop over all the satellite pixels and identify which GCAS pixels lie
% within each one. We'll just check if the center of the GCAS pixel is
% inside the OMI pixel because the GCAS ones are significantly smaller

for a=1:numel(this_data.Longitude)
    xx = inpolygon(GCAS.Longitude, GCAS.Latitude, this_data.FoV75CornerLongitude(:,a), this_data.FoV75CornerLatitude(:,a));
    this_data.Matches(a) = sum(xx) > 0;
    if this_data.Matches(a)
        for f=1:numel(fns)
            new_fn = sprintf('GCAS_%s',fns{f});
            this_data.(new_fn)(a) = nanmean(GCAS.(fns{f})(xx));
        end
    end
end

end
