function [  ] = load_all_merges( merge_dir, year_flag )
%load_all_merges Loads all merge .mat files from a specified directory.
%   Loads all merge files from a directory specified, either as a string
%   input to this function, or via an open dialogue presented to the user.
%
%   Each merge file will be loaded as a structure with the variable name
%   "MergeMMDD" where MMDD is the numeric month and day the file comes
%   from. Set the second input to true (or a non-zero scalar) in order to
%   make the variable names MergeYYYYMMDD - if you will be loading multiple
%   campaigns and need to keep them straight. If you wish to have the names
%   printed with year but still use the UI to choose the folder to load,
%   simply pass a boolean or scalar as the only argument - as long as it is
%   not a string, this function will recognize your intent.

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INITIALIZATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

E = JLLErrors;
narginchk(0,2);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    year_flag = false;
end


if nargin == 1 && ~ischar(merge_dir)
    % If the user does not pass a directory to load, open a UI to load the
    % folder from. But, if something is passed other than a string, assume
    % that it is the year_flag variable and thus set it, overriding the
    % assumed false value from above.
    
    if ~ischar(merge_dir)
        if nargin > 1;
            E.badinput('If you are passing two inputs, the first must be a string specifying the directory of the merge files and the second a boolean specifying whether to include the year in the var names');
        else % If the user only passes one input assume that it is the year_flag, not the merge directory.
            year_flag = merge_dir;
            merge_dir = uigetdir('/Volumes','Select the directory with the Merge files in it');
        end
    end
elseif nargin == 0
    merge_dir = uigetdir('/Volumes','Select the directory with the Merge files in it');
end

%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN BODY %%%%%
%%%%%%%%%%%%%%%%%%%%%

% Identify all .mat files in the merge directory
files = dir(fullfile(merge_dir,'*.mat'));
for i=1:numel(files)
    % Load each file
    filename = fullfile(merge_dir, files(i).name);
    M = load(filename, 'Merge');
    
    % Now find the date, which should be saved in the metadata of the Merge
    % structure. If the date field doesn't exist, error out. Use datestr to
    % ensure formatting as YYYY-MM-DD
    if ~isfield(M.Merge,'metadata') || ~isfield(M.Merge.metadata,'date')
        E.callError('date_not_found','Could not identify the date field in the Merge structure - there should be Merge.metadata.date in the file to be opened');
    end
    filedate = datestr(M.Merge.metadata.date,29);
    if year_flag
        varname = sprintf('Merge%s%s%s',filedate(1:4),filedate(6:7),filedate(9:10));
    else
        varname = sprintf('Merge%s%s',filedate(6:7),filedate(9:10));
    end
    
    % Save out the variable in the base workspace, then clear (no sense in
    % wasting memory)
    eval(sprintf('%s = M.Merge',varname));
    putvar(sprintf('%s',varname));
    clear(varname);
end

end

