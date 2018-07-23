function [ behr_flags, flags_meaning ] = behr_quality_flags( data )
%BEHR_QUALITY_FLAGS Create the BEHRQualityFlags field
%   In order to simplify things for the end user, I decided to combine the
%   NASA quality flags with some of our own into a single flag field. Like
%   VcdQualityFlags and XTrackFlags, the idea is that we will use a bit
%   array, i.e. a number where each bit that makes up its binary
%   representation has a specific meaning.
%
%   The first (least significant) bit is reserved as an QUALITY summary
%   flag, i.e. it will be == 1 if a problem with the pixel means that it
%   should not be used when trying to retrieve a to-ground (visible+ghost)
%   tropospheric VCD. The idea is that end users will be able to remove any
%   pixel where this bit is set (which means the flag value will be odd)
%   and that will be sufficient filtering to remove any low quality pixels.
%
%   The second bit is reserved as an ERROR summary flag. It will be set to
%   1 if any more specific error flag is set. This indicates that the pixel
%   should not be used under ANY case for a to-ground or visible-only
%   tropospheric VCD.
%
%   The third through sixteenth bits are to provide more information about
%   the root cause of the error, each one should be set if a certain
%   condition is met. If any of these are set, then the first bit (summary
%   bit) should also be set.
%
%   The seventeenth through thirty second bits indicate warnings, these
%   could be something that indicates that there MIGHT be a problem with
%   the pixel, but that it is usually still useable, or just noting that a
%   certain behavior occured in the retrieval for that pixel; for instance,
%   I intend to use this to mark instances where the BRDF albedo used a
%   water model, rather than the MODIS parameters.
%
%   Usage:
%
%   BEHR_FLAGS = BEHR_QUALITY_FLAGS( DATA ) where DATA is a scalar
%   structure containing the native resolution BEHR data. This function
%   requires it have the fields:
%       * BEHRAMFTrop
%       * BEHRAMFTropVisOnly
%       * VcdQualityFlags
%       * XTrackQualityFlags
%       * AlbedoOceanFlag
%       * MODISAlbedoQuality
%       * MODISAlbedoFillFlag
%       * CloudFraction
%       * TropoPresVSCldPres
%       * InterpolatedTropopausePressure
%   [ ___, FLAG_MEANING ] = BEHR_QUALITY_FLAGS( DATA )
%   [  ~ , FLAG_MEANING ] = BEHR_QUALITY_FLAGS() The second output is a
%   cell array giving the meaning of each flag bit. This is output either
%   if DATA is given as an input or if no inputs are given. In the latter
%   case, the value for the BEHR_FLAGS output should be discarded, since it
%   is just a fill value since no data to base the flags on was given.

E = JLLErrors;


% Add any fields used to set the flags here. They will be checked that they
% are in data, or used to set default values if the user just needs the
% flag meaning.
req_fields = {'BEHRAMFTrop', 'BEHRAMFTropVisOnly', 'VcdQualityFlags', 'XTrackQualityFlags',...
    'AlbedoOceanFlag', 'MODISAlbedoQuality', 'MODISAlbedoFillFlag', 'CloudFraction', 'TropoPresVSCldPres','Interp_TropopausePressure'};

if nargin == 0    
    % If given no arguments, we must just want the flag meanings, so create
    % filler arrays to allow the set_flags calls to work. Assign a logical
    % type b/c some of the set_flags read the fields directly and expect
    % that input to be a logical array.
    for a=1:numel(req_fields)
        data.(req_fields{a}) = false;
    end
else
    if ~isstruct(data) || ~isscalar(data)
        E.badinput('DATA must be a scalar structure');
    end

    xx = ~ismember(req_fields, fieldnames(data));
    if any(xx)
        E.badinput('DATA is missing the required field(s): %s', strjoin(req_fields(xx), ', '));
    end
end

% Set up the output arrays and define the summary bits

% these are shared with the nested function set_flags, so do not need to be
% passed to or returned from that function. This is a rare case where it
% makes sense to not pass the variable explicitly, because "set_flags"
% clearly should modify the flags.
behr_flags = uint32(zeros(size(data.BEHRAMFTrop)));
flags_meaning = cell(1, 32);
flags_meaning{1} = '1: Summary bit: set if pixel would produce poor quality to-ground VCD';
flags_meaning{2} = '2: Critical error bit: set if VCD from pixel should not be used under any condition';

%%%%%%%%%%%%%%%
% ERROR FLAGS %
%%%%%%%%%%%%%%%

% Set an error flag if the AMF has been set to the minimum value
set_flags(data.BEHRAMFTrop <= behr_min_amf_val() | isnan(data.BEHRAMFTrop)...
    | data.BEHRAMFTropVisOnly <= behr_min_amf_val() | isnan(data.BEHRAMFTropVisOnly),...
    3, true, true, 'BEHR AMF error: AMF below minimum value');

% Set an error flag if the VcdQualityFlags field is not an even value (it's
% own quality summary flag was set)
set_flags(mod(data.VcdQualityFlags, 2) ~= 0, 4, true, true, 'VcdQualityFlags: NASA summary flag set');

% Set an error flag if the XTrackQualityFlags fields is > 0, i.e. it has
% been affected by the row anomaly. Had to switch away from ~= 0 because
% NaN ~= 0 will return true, but NaN > 0 will return false. For early files
% (~2006 to ~Nov 2008 last I checked), the XTrackQualityFlags fields seems
% to have all fill values, for orbits not affected by the row anomaly, so
% this relies on fill values being imported as NaNs and thus not triggering
% this logical test.
set_flags(data.XTrackQualityFlags > 0, 5, true, true, 'XTrackQualityFlags: NASA flag > 0');


%%%%%%%%%%%%%%%%%
% WARNING FLAGS %
%%%%%%%%%%%%%%%%%

% Set a warning flag if the OMI cloud fraction is > 20%. This will also
% indicate that the pixel should not be used for to-ground VCDs.
set_flags(data.CloudFraction > 0.2, 17, true, false, 'OMI effective geometric cloud fraction >20%');

% Set a warning flag if we have to use an ocean model for the BRDF
set_flags(data.AlbedoOceanFlag, 18, false, false, 'Ocean Albedo Flag: surface albedo uses COART LUT');

% Set a warning if the BRDF quality is worse that 2.5. A quality of 2 is
% "relative good quality, 75% or more with full inversions" while 3 is
% "mixed, <= 75% full inversions and <= 25% fill values". Having the cut
% off be 2.5 will all some poor quality MODIS BRDFs to contribute to the
% surface reflectance without automatically flagging it as poor quality.
set_flags(data.MODISAlbedoQuality >= 2.5 | data.MODISAlbedoFillFlag, 19, true, false, 'MODIS BRDF quality worse than ( >= ) 2.5 or > 50% of MODIS grid cells had fill value');

% Set a warning if the cloud pressure is smaller than tropopause pressure
set_flags(data.TropoPresVSCldPres, 20, false, false, 'Cloud pressure is smaller than tropopause pressure');

% Set a warning if the wrf tropopause pressure is not calculated by lapse
% rate but interpolated from neighboring points
set_flags(data.Interp_TropopausePressure, 21, false, false, 'Tropopause Presssure is interpolated from neighboring points as no point with a lapse rate < 2 K/km is found');

% NOTE: the 32nd bit is reserved for fill values. This is defined in
% BEHR_publishing_attribute_table and used in BEHR_publishing_main when
% creating the HDF files. Although these quality flags should never have a
% fill value, it is good practice to leave one available. Therefore, the
% 32nd bit should NEVER be used here.

    function set_flags(bool_mask, bit, bad_to_ground_quality, is_error, explanation_string)
        % This nested subfunction should always be used to set the flags.
        %
        %   BOOL_MASK: a logical array the same size as FLAGS that is true
        %   where the bit should be set in the FLAGS array.
        %
        %   BIT: the index of the bit to set to 1 (1-based, starting from
        %   least significant)
        %
        %   BAD_TO_GROUND_QUALITY: a scalar logical that, if true,
        %   indicates that this bit is a reason that the pixel would give
        %   bad quality to-ground tropospheric VCDs. It will cause the
        %   summary bit to be set everywhere that this is true. Unlike
        %   WARNING_ONLY (below), this has no relation to which bit
        %   position is being modified, since it is entirely possible that
        %   a warning bit, which doesn't indicate a serious flaw in the
        %   retrieval, could be used to indicate that, for whatever reason,
        %   this particular pixel wouldn't produce good quality to-ground
        %   (so including ghost column) tropospheric VCDs. An example is
        %   that cloud fraction exceeds some threshold.
        %
        %   IS_ERROR: a scalar logical that, if true, indicates that this
        %   bit is an wrror flag. When true, it will set the critical error
        %   and summary bits to 1 everywhere that BIT is also set to true.
        %   It is intended that bits 3-16 are error bits (i.e. IS_ERROR =
        %   true) and 17-32 are warning bits (IS_ERROR = false). A warning
        %   will be issued if the bit number and IS_ERROR do not agree,
        %   this is a check to make sure any future bits follow this
        %   convention.

        % Type and size checking
        if ~isinteger(behr_flags)
            E.badinput('FLAGS must be an integer type')
        elseif ~islogical(bool_mask) || ~isequal(size(bool_mask), size(behr_flags))
            E.badinput('BOOL_MASK must be a logicial type the same size as FLAGS')
        elseif ~isscalar(bit) || ~isnumeric(bit) || bit < 1
            E.badinput('BIT must be a positive scalar number')
        elseif ~isscalar(bad_to_ground_quality) || ~islogical(bad_to_ground_quality)
            E.badinput('BAD_TO_GROUND_QUALITY must be a scalar logical');
        elseif ~isscalar(is_error) || ~islogical(is_error)
            E.badinput('IS_ERROR must be a scalar logical');
        end
        
        if bit < 3
            E.badinput('BIT < 3: the first two bits are reserved for summary flags')
        elseif (bit > 16 && is_error) || (bit <= 16 && ~is_error)
            warning('It is expected that bits 2-16 are error bits and 17-32 are warning bits; if you are following this, be sure WARNING_ONLY reflects it');
        end
        
        behr_flags = bitset(behr_flags, bit, bool_mask);
        if is_error
            error_mask = bool_mask | bitand(behr_flags, 2);
            behr_flags = bitset(behr_flags, 2, error_mask);
        end
        if is_error || bad_to_ground_quality
            bad_qual_mask = bool_mask | bitand(behr_flags,1);
            behr_flags = bitset(behr_flags, 1, bad_qual_mask);
        end
        flags_meaning{bit} = sprintf('%d: %s', bit, explanation_string);
    end

end
