function [  ] = examine_aerosol_categorization( campaign_name, profaltmat, varargin )
%examine_aerosol_categorization Plots the LIF and NCAR NO2 profiles against
%aerosol extinction profiles with the 90% column heights marked. Takes 2
%arguments: a campaign name as a string and one of the altitude fields of
%the structure output from categorize_aerosol_profiles.
%
%   Parameters:
%       norm = Defaults to true. Boolean whether to normalize profiles to
%       their maximum value.

DEBUG_LEVEL = 1;

% Input checking and parsing
E = JLLErrors;
narginchk(2,2);
if ~ischar(campaign_name);
    error(E.badinput('''campaign_name'' must be a string.'));
elseif isstruct(profaltmat)
    error(E.badinput('''profaltmat'' should be a matrix - one of the altitude fields from the categorization structure.  Did you pass the whole structure?'));
elseif ~isnumeric(profaltmat)
    error(E.badinput('''profaltmat'' should be a numeric matrix'))
end

% set the default values of the parameters
norm_bool = 1;

parameters = varargin;
while numel(parameters) > 0
    par_name = lower(parameters{1});
    par_val = parameters{2};
    if ~ischar(par_name)
        error(E.badinput,'Parameter names must be strings');
    end
    switch par_name
        case 'norm'
            if ~isscalar(par_val)
                error(E.badinput('The parameter ''norm'' must be a scalar value'));
            end
            norm_bool = par_val;
        otherwise
            error(E.badinput(sprintf('"%s" is not a recognized parameter',par_name)));
    end
    parameters(1:2) = [];     
end

% Load the campaign data based on the campaign name
[ Names, date_cell, directory, range_file ] = merge_field_names( campaign_name );

% Set the start and end times you'll consider. use military time
starttime = '10:30';
endtime = '16:30';
% The timezone.  Set to 'auto' to calculate from longitude.
tz = 'auto';



profilefield = Names.profile_numbers;
no2field = Names.no2_lif;
no2ncarfield = Names.no2_ncar;
aerosol_field = Names.aerosol_extinction;
alt_field = Names.gps_alt;

% This pattern is used to load the day's merge file; it should not need
% edited to be more specific than '*%s_%s_%s.mat'.
mergepat = '*%s_%s_%s.mat';

% Load the range file if there are no profile numbers. We will use the
% boolean set here to control how the program IDs profiles - by profile
% number or UTC range.
if ~isempty(profilefield)
    prof_bool = true;
else
    SRange = load(range_file);
    RangesStruct = SRange.Ranges;
    range_dates = cellstr(datestr({RangesStruct.Date},29));
    prof_bool = false;
end

if strcmpi(tz,'auto');
    tz_bool = true;
else
    tz_bool = false;
end

dates = datenum(date_cell{1}):datenum(date_cell{2});
for d=1:numel(dates)
    curr_date = datestr(dates(d),29);
    if DEBUG_LEVEL > 0; fprintf('Now checking %s\n',curr_date); end
    
    % Load the merge file and extract the necessary data
    SNO2 = wildcard_load(directory,mergepat,dates(d),'Merge');
    if isempty(SNO2); continue; end
    Merge = SNO2.Merge;
    
    [no2, utc, alt, lon] = remove_merge_fills(Merge,no2field,'alt',alt_field);
    no2ncar = remove_merge_fills(Merge, no2ncarfield);
    ext = remove_merge_fills(Merge,aerosol_field);
    
    % If the timezone was set to "auto," calculate the difference from UTC
    % based on the mean longitude. This will produce a vector of the same
    % length as the utc and lon variables
    if tz_bool
        tz = round(lon/15);
    end
    if prof_bool
        prof_bool = true;
        profnum = remove_merge_fills(Merge,profilefield);
        % Identify the profiles that fall within the start and end times. 
        
        % Now, use the timezone (entered or calculated) to produce a new vector of
        % times that reflect the local time at each point
        local_times = utc2local_sec(utc,tz);
        
        % Get all unique profile numbers and their start times
        unique_profnums = unique(profnum(profnum~=0));
        start_times = zeros(numel(unique_profnums),1);
        for a=1:numel(unique_profnums)
            % If we're using a vector of timezones, get the most common
            % timezone from the profile - we'll treat the whole profile as
            % belonging to that timezone.  If timezone has been set
            % manually, then just use that time zone to convert the start
            % time of the profile. (obviously) Save the local start time.
            xx = profnum == unique_profnums(a);
            if ismatrix(tz) && isnumeric(tz)
                % Case where we're using a vector of timezones
                mct = mode(tz(xx));
                start_times(a) = utc2local_sec(min(utc(xx)),mct);
            elseif ischar(tz)
                start_times(a) = utc2local_sec(min(utc(xx)),tz);
            else
                error(E.callError('tz_not_recognized','Cannot understand the format of timezones'));
            end
        end
        
        % Remove from consideration any profiles with a start time before 10:45
        % am or after 4:45 pm local standard time
        yy = start_times >= local2utc(starttime,0) & start_times <= local2utc(endtime,0);
        unique_profnums = unique_profnums(yy); start_times = start_times(yy);
        num_profs = numel(unique_profnums);
    else
        % Find the set of ranges that correspond to this date
        rr = find(strcmp(curr_date,range_dates));
        if isempty(rr);
            error(E.callError('ranges',sprintf('No UTC ranges found for %s',curr_date)));
        end
        Ranges = RangesStruct(rr).Ranges;
        
        % Find the ranges that fall in the start and end times
        yy = false(size(Ranges,1),1);
        for a=1:size(Ranges,1)
            tz_ind = utc >= Ranges(a,1) & utc <= Ranges(a,2);
            mct = mode(tz(tz_ind));
            range_start_local = utc2local_sec(Ranges(a,1),mct);
            yy(a) = range_start_local >= local2utc(starttime,0) && range_start_local <= local2utc(endtime,0);
        end
        
        ranges_in_time = Ranges(yy,:);
        num_profs = size(ranges_in_time,1);
    end
    
    for pn = 1:num_profs;
        % We need to create the logical index matrix pp differently whether
        % we're using profile numbers or UTC ranges
        if prof_bool
            p = unique_profnums(pn);
            % Skip this loop if none of the profiles are in the structure
            % output from categorize_aerosol_profiles.
            if all(p~=profaltmat(:,1)); continue; end
            
            % Find which altitudes correspond to the profiles
            ff = p == profaltmat(:,1);
            % Save the altitude values; which column to use is different
            % for NO2 and aerosols, and whether the IDs are profile numbers
            % or UTC ranges.
            no2_90alts = profaltmat(:,2);
            aer_90alts = profaltmat(:,3);
            
            % Finally set the logical index matrix
            pp = profnum == p;
        else
            p = ranges_in_time(pn,:);
            
            % Skip this if the range isn't in the indicated field. Do this
            % by only comparing the start times.
            if all(p(1,1) ~= profaltmat(:,1)); continue; end
            
            % Find which altitudes correspond to the profiles
            ff = p(1,1) == profaltmat(:,1) & p(1,2) == profaltmat(:,2);
            % Save the altitude values; which column to use is different
            % for NO2 and aerosols, and whether the IDs are profile numbers
            % or UTC ranges.
            no2_90alts = profaltmat(:,3);
            aer_90alts = profaltmat(:,4);
            
            % Finally set the logical index matrix
            pp = utc >= p(1,1) & utc <= p(1,2);
        end
        
        % First bin both LIF and NCAR NO2 data plus the aerosol extinction
        % data.
        [no2bins, no2alts] = bin_rolling_vertical_profile(alt(pp),no2(pp),0.5,0.1);
        [ncarbins, ncaralts] = bin_rolling_vertical_profile(alt(pp), no2ncar(pp), 0.5, 0.1);
        [aerbins, airalts] = bin_rolling_vertical_profile(alt(pp),ext(pp),0.5,0.1);
        
        % Axis labels with units. Overwritten if bins are normalized.
        ax_label = {'Aerosol extinction at 532 nm','NO2 pptv'};
        
        % Normalize the profiles to their maximum value if indicated to do
        % so by the user (true by default)
        if norm_bool
            no2bins = no2bins ./ max(no2bins);
            ncarbins = ncarbins ./ max(ncarbins);
            aerbins = aerbins ./ max(aerbins);
            
            ax_label = {'Aerosol extinction at 532 nm (normalized)','NO2 (normalized)'};
        end
        
        % Plot the LIF and aerosol data initially.
        figure;
        hax = plotxx(aerbins, airalts, no2bins, no2alts, ax_label,{'Alt (km)','Alt (km)'});
        % Then reset the limits to include 0 if the aerosol data doesn't
        % cross 0. Plot the aerosol 90% height.
        oldlim = get(hax(1),'xlim'); 
        newlim = [min(0,min(aerbins)),max(0,max(oldlim))];
        set(hax(1),'xlim', newlim);
        line([0, oldlim(2)/2, oldlim(2)],[aer_90alts(ff),aer_90alts(ff),aer_90alts(ff)],'linestyle','--','color','k','linewidth',2,'marker','x','parent',hax(1)); 
        
        % Now again set the limits to include 0 but for the NO2 axes.
        oldlim = get(hax(2),'xlim'); 
        newlim = [min(0,min(no2bins)),max(0,max(oldlim))];
        set(hax(2), 'xlim', newlim);
        % Add the NCAR data
        line(ncarbins, ncaralts, 'color', [0 0.7 0], 'parent',hax(2));
        % And draw the NO2 90% altitude.
        line([0, oldlim(2)],[no2_90alts(ff),no2_90alts(ff)],'linestyle','--','color','r','linewidth',2,'parent',hax(2));
        
        ylims1 = get(hax(1),'ylim');
        ylims2 = get(hax(2),'ylim');
        newylim = [0, max([ylims1, ylims2])];
        set(hax(1),'ylim',newylim);
        set(hax(2),'ylim',newylim);
        
        dt = Merge.metadata.date;
        title(sprintf('%s Prof Num %d',dt,p));
    end
    if ~isempty(get(0,'children'))
        tilefigs;
        pause;
        close all
    end
end


end

