function [  ] = resave_mats(  )
%resave_mats Loads and resaves .mat files in a folder given.
%   I've been having trouble where .mat files sent to me will cause a new
%   instance of MATLAB to open whenever I double click on one outside of
%   MATLAB.  The fix seems to be to load them once and resave a new mat
%   file myself.  This may be related to Mac permission, but I've no idea
%   why it started when I updated to Matlab 2014b from 2013a.

THE_MOST_SPECIAL_FOLDER_EVER = uigetdir;

ALL_MY_FILES_DAVE = dir(fullfile(THE_MOST_SPECIAL_FOLDER_EVER,'*.mat'));

% We're going to be loading variables from the files, so we need ALL of the
% non-loaded variables present before we start loading things so that we
% know which variables to save and clear.
%
% I'm deliberately using absurd variable names to minimize the chances of a
% loaded variable matching one of the variables in this function.
CLEAR_ALL_THE_THINGS = {};
CRAZY_INDEX_PATRICE=1;
ORIGINAL_VARS_allen = struct;
NEW_VARS_beatrice = struct;
ORIGINAL_VARS_allen = whos;

for CRAZY_INDEX_PATRICE=1:numel(ALL_MY_FILES_DAVE);
    load(fullfile(THE_MOST_SPECIAL_FOLDER_EVER,ALL_MY_FILES_DAVE(CRAZY_INDEX_PATRICE).name));
    NEW_VARS_beatrice = whos;
    CLEAR_ALL_THE_THINGS = find_new_vars(ORIGINAL_VARS_allen,NEW_VARS_beatrice);
    save(fullfile(THE_MOST_SPECIAL_FOLDER_EVER,ALL_MY_FILES_DAVE(CRAZY_INDEX_PATRICE).name),CLEAR_ALL_THE_THINGS{:});
    clear(CLEAR_ALL_THE_THINGS{:});
end
end

function newvars = find_new_vars(W1,W2)
    vars1 = {W1(:).name};
    vars2 = {W2(:).name};
    xx = ~containedin(vars2,vars1);
    newvars = vars2(xx);
end