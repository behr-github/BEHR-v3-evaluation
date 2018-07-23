function [ F_behr, behr_dir ] = list_behr_files( start_date, end_date, varargin )
%LIST_BEHR_FILES List all BEHR .mat files between two dates
%   [ F, BEHR_DIR ] = LIST_BEHR_FILES( START_DATE, END_DATE ) will list all
%   monthly profile, US BEHR files between START_DATE and END_DATE, which
%   must be date specifications understood by Matlab (strings or numbers).
%   Returns F, a structure from DIR(), and BEHR_DIR, the directory it
%   searched.
%
%   [ ___ ] = LIST_BEHR_FILES( ___, PROF_MODE )
%   [ ___ ] = LIST_BEHR_FILES( ___, PROF_MODE, REGION ) 
%   [ ___ ] = LIST_BEHR_FILES( ___, PROF_MODE, REGION, VERSION ) allows you
%   to use different profile mode BEHR files ('daily' or 'monthly'),
%   different regions ('us', or 'hk'), and different versions.
%
%   [ ___ ] = LIST_BEHR_FILES( ___, 'all' ) will list all BEHR files for
%   the given date range(s), even if they do not exist.

E = JLLErrors;

p = advInputParser;
p.addOptional('prof_mode', 'monthly', @ischar, 'PROF_MODE must be a string');
p.addOptional('region', 'us', @ischar, 'REGION must be a string');
p.addOptional('version', BEHR_version(), @ischar, 'VERSION must be a string');
p.addFlag('all');

p.parse(varargin{:});
pout = p.Results;

prof_mode = pout.prof_mode;
region = pout.region;
version_str = pout.version;
all_files = pout.all;

start_date = validate_date(start_date);
end_date = validate_date(end_date);

if numel(start_date) ~= numel(end_date)
    E.badinput('Numel of elements in START_DATE and END_DATE must be equal');
end

F_behr = [];
file_pattern = behr_filename('*', prof_mode, region, '.mat', version_str);

for a=1:numel(start_date)
    if a == 1
        behr_dir = behr_paths.BEHRMatSubdir(region, prof_mode);
    else
        if ~strcmp(behr_dir, behr_paths.BEHRMatSubdir(region, prof_mode))
            E.notimplemented('LIST_BEHR_FILES cannot handle BEHR files across multiple directories');
        end
    end
    F = dir(fullfile(behr_dir, file_pattern));
    
    if isempty(F)
        continue
    end
    
    dvec = datenum(regexp({F.name}, '\d\d\d\d\d\d\d\d', 'match', 'once'), 'yyyymmdd');
    dd = dvec >= start_date(a) & dvec <= end_date(a);
    
    F(~dd) = [];
    if all_files
        F = insert_missing_dates(F, start_date(a), end_date(a));
    end
    
    F_behr = veccat(F_behr, F);

end

    function F_full = insert_missing_dates(F_in, start_date_in, end_date_in)
        full_dvec = start_date_in:end_date_in;
        files_dvec = date_from_behr_filenames(F_in);
        xx_present = ismember(full_dvec, files_dvec);
        F_full(xx_present) = F_in;
        missing_inds = find(~xx_present);
        for i_missing = missing_inds(:)' % must be row vector so that i_missing takes on each valing in missing_inds
            F_full(i_missing).name = behr_filename(full_dvec(i_missing), prof_mode, region, '.mat', version_str);
        end
    end

end

