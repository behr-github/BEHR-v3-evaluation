function [file, filter_idx] = ask_file(prompt, file_or_dir, varargin)
% ASK_FILE Prompt the user to specify a file or directory
%   [ FILE, FILTER_IDX ] = ASK_FILE( PROMPT, FILE_OR_DIR ) Will prompt the
%   user to give a file or directory. PROMPT must be a string that will be
%   printed to prompt the user to give the file or directory. FILE_OR_DIR
%   must be one of the strings 'file', 'dir', or 'either', controlling
%   whether the user is expected to provide a file, directory, or either.
%
%   How this function behaves depends on two things: whether the
%   MATLAB_DISPLAY environmental variable is set to 0 and whether
%   FILE_OR_DIR is 'either'. If MATLAB_DISPLAY is 0 or FILE_OR_DIR is
%   'either' then the user will be prompted to enter a path as text.
%   Otherwise, a file dialogue is opened using UIGETFILE or UIGETDIR. (If
%   FILE_OR_DIR is 'either', we have to do things by text because UIGETFILE
%   cannot select directories and UIGETDIR cannot selet files.)
%
%   The output FILE will be the full path to the file or directory selected
%   and, if FILE_OR_DIR is 'file' and the UI selector is used then
%   FILTER_IDX will be set to the return value FILTERINDEX from UIGETFILE.
%   Otherwise it will be a NaN.
%
%   Parameters:
%
%       'softquit' - if false (default) a user cancelled error is thrown if
%       the user aborts selecting a file. If true, then FILE will be
%       returned as 0 in that case.
%
%       'force_cl' - if true, will always use the command line input, even
%       if the UI could be used. Default is false.
%
%   Other optional and parameter args will be passed on to UIGETFILE or
%   UIGETDIR if used. Note that optional arguments (e.g. FILTERSPEC,
%   TITLE, or FILE for UIGETFILE) must come before any parameter args,
%   whether the parameters are for ASK_FILE or UIGETFILE/UIGETDIR.

E = JLLErrors;
p = inputParser;

% These are actual parameters for this function
p.addParameter('softquit', false);
p.addParameter('force_cl', false);
p.addParameter('ui_args', {});

p.parse(varargin{:});
pout = p.Results;

soft_quit = pout.softquit;
force_cl = pout.force_cl;
ui_args = pout.ui_args;

allowed_file_or_dir = {'file','dir','either'};

if ~ischar(prompt)
    E.badinput('PROMPT must be a char array');
elseif ~ischar(file_or_dir) || ~ismember(file_or_dir, allowed_file_or_dir)
    E.badinput('FILE_OR_DIR must be one of the following char arrays: %s', strjoin(allowed_file_or_dir, ', '));
end

if ~islogical(soft_quit) || ~isscalar(soft_quit)
    E.badinput('''soft_quit'' must be a scalar logical');
elseif ~islogical(force_cl) || ~isscalar(force_cl)
    E.badinput('''force_cl'' must be a scalar logical');
end

if ~force_cl && ~strcmpi(file_or_dir, 'either') && isDisplay
    % If no filterspec is given, force it to default to all files, if we're
    % selecting a file. Otherwise the UI allows only a limited selection of
    % files
    if strcmpi(file_or_dir, 'file') && isempty(ui_args) 
        ui_args = veccat({'*'}, ui_args);
    end
    [file, filter_idx] = use_ui(prompt, file_or_dir, ui_args);
else
    file = use_cl(prompt, file_or_dir);
    filter_idx = nan;
end

if ~soft_quit && isnumeric(file) && file == 0
    E.userCancel();
end

end


function [ file, filter_idx ] = use_ui(prompt, file_or_dir, ui_args)
% print the prompt to the command window before bringing up the UI because
% sometimes the UI window doesn't show a very useful prompt, even if given
% a title (I've had this issue on Macs)
fprintf('%s\n', prompt);
input('Press ENTER to continue');

if strcmpi(file_or_dir, 'file')
    [filenames, pathname, filter_idx] = uigetfile(ui_args{:});
    if isnumeric(filenames) % user cancelled, filenames should be 0
        file = filenames;
    elseif iscell(filenames)
        file = cell(size(filenames));
        for i = 1:numel(file)
            file{i} = fullfile(pathname, filenames{i});
        end
    else
        file = fullfile(pathname, filenames);
    end
else
    file = uigetdir(ui_args{:});
    filter_idx = nan;
end

end

function file = use_cl(prompt, file_or_dir)
E = JLLErrors;

switch lower(file_or_dir)
    case 'file'
        % exist(f, 'file') will return 7 if f is a directory, 2 if it is a
        % proper file
        check_fxn = @(f) exist(f, 'file') == 2;
        error_str = 'file';
    case 'dir'
        check_fxn = @(f) exist(f, 'dir');
        error_str = 'directory';
    case 'either'
        check_fxn = @(f) exist(f, 'file');
        error_str = 'file or directory';
    otherwise
        E.notimplemented('No method for file_or_dir == "%s"', file_or_dir);
end

while true
    file = input(sprintf('%s: ', prompt), 's');
    if strcmpi(file, 'q')
        file = 0;
        return
    elseif check_fxn(file)
        return
    else
        if isempty(file)
            file_error = 'An empty string';
        else
            file_error = file;
        end
        fprintf('%s is not a valid %s. Please try again, or type "q" to quit', file_error, error_str);
    end
end
end
