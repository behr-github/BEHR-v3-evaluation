function behr_dates = date_from_behr_filenames(behr_files)
%DATE_FROM_BEHR_FILENAMES Generate a list of date numbers from BEHR files
%   BEHR_DATES = DATE_FROM_BEHR_FILENAMES( BEHR_FILES ) Given a single BEHR
%   files, a cell array of BEHR file names, or a structure of files
%   returned by DIR(), generates a vector of date numbers corresponding to
%   the date in each of the file names.

behr_files = files_input(behr_files);
behr_dates = nan(size(behr_files));
for i=1:numel(behr_files)
    behr_dates(i) = datenum(regexp(behr_files{i}, '\d{8}(?=\.mat|\.hdf)', 'match', 'once'), 'yyyymmdd');
end

end

