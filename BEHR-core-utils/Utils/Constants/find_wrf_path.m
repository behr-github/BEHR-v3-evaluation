function [ wrf_path ] = find_wrf_path( region, profile_mode, this_date, varargin )
%FIND_WRF_PATH Returns the path to WRF profiles for the given date
%   WRF_PATH = FIND_WRF_PATH( REGION, PROFILE_MODE, THIS_DATE ) Returns the
%   proper path to the WRF-Chem profiles as WRF_PATH. Given PROFILE_MODE =
%   'monthly', will just return behr_paths.wrf_monthly_profiles. Given
%   PROFILE_MODE = 'daily', it will search all paths defined in the cell
%   array behr_paths.wrf_profiles for the year and month subdirectories
%   that match THIS_DATE. If it finds one, it will return it. Otherwise, an
%   error is thrown. It does not verify that the required wrfout files are
%   actually present. THIS_DATE must be either a date number or date string
%   implicitly understood by Matlab.
%
%   WRF_PATH = FIND_WRF_PATH( REGION, PROFILE_MODE, THIS_DATE, 'fullpath' )
%   will include the file corresponding to the given date. Note that it
%   offers no guarantee that file exists; it will return the file for the
%   exact hour/minute/second requested.

E = JLLErrors;

p = advInputParser;
p.addFlag('fullpath');

p.parse(varargin{:});
pout = p.AdvResults;

find_exact_file = pout.fullpath;

this_date = validate_date(this_date);

if strcmpi(profile_mode, 'monthly')
    wrf_path = behr_paths.wrf_monthly_profiles;
    if find_exact_file
        wrf_path = fullfile(wrf_path, monthly_file_name(this_date));
    end
    return
elseif strcmpi(profile_mode, 'daily')
    yr_str = datestr(this_date, 'yyyy');
    mn_str = datestr(this_date, 'mm');
    wrf_dirs = behr_paths.wrf_profiles;
    for a=1:numel(wrf_dirs)
        if ~exist(wrf_dirs{a}, 'dir')
            E.dir_dne('The root WRF-Chem profile directory %s does not exist;\n  Is the file server mounted?\n  Is the directory defined correctly in behr_paths.m?', wrf_dirs{a});
        end
        
        wrf_path = fullfile(wrf_dirs{a}, region, yr_str, mn_str);
        if exist(wrf_path, 'dir')
            if find_exact_file
                wrf_path = fullfile(wrf_path, daily_file_name(this_date));
            end
            
            return
        end
    end
    
    % If we've gotten here, we haven't found the directory
    E.dir_dne('No WRF-Chem daily output directory exists for %s in region %s', datestr(this_date, 'mmm yyyy'), upper(region));
else
    E.badinput('No paths defined for PROFILE_MODE = ''%s''', profile_mode);
end

end

function filename = monthly_file_name(this_date)
filename = sprintf('WRF_BEHR_monthly_%s.nc', datestr(this_date, 'mm'));
end

function filename = daily_file_name(this_date)
% Need to round this_date to the nearest hour
time_component = mod(this_date,1)*24;
this_date = floor(this_date) + round(time_component)/24;
filename = sprintf('wrfout_d01_%s', datestr(this_date, 'yyyy-mm-dd_HH-00-00'));
end