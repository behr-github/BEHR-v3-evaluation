%record_profile_ranges
%
%   This script will run through all the merge files in the given directory
%   and call the "select_changing_altitude" GUI to allow the user to mark
%   off the UTC ranges that correspond to parts of the aircraft's flight
%   where substantial altitude changes occur that will provide a complete
%   vertical profile.  The resulting UTC ranges will be stored in a data
%   structure called "Ranges" where each top level index corresponds to a
%   file; the date and filename will be recorded along with the table of
%   UTC ranges.

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%     INPUT    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

% This is the directory where the merge files to be analyzed are.  Leave it
% as an empty string to select the directory with a dialog
merge_dir = '';

if isempty(merge_dir);
    merge_dir = uigetdir('/Volumes','Select the directory with the merge files in it');
end

% Add the path with the GUI in it
addpath('/Users/Josh/Documents/MATLAB/NO2 Profiles/GUIs');

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%    MAIN LOOP    %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find all the matlab files in the directory; iterate through them, loading
% the file into the GUI, then save the information 
merge_files = dir(fullfile(merge_dir,'*.mat'));
n = numel(merge_files);
Ranges = struct('Date',cell(1,n),'Filename',cell(1,n),'Ranges',cell(1,n));
for a=1:n
    filename = fullfile(merge_dir, merge_files(a).name);
    r = select_changing_altitude(filename);
    load(filename,'Merge');
    Ranges(a).Date = Merge.metadata.date;
    Ranges(a).Filename = filename;
    Ranges(a).Ranges = r;
end