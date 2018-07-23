function [ repo_dir ] = no2_prof_repo(  )
%NO2_PROF_REPO Return the absolute path to the NO2 Profile repository

mydir = fileparts(mfilename('fullpath'));
%Assumes this is in the utility_scripts directory within the BEHR repo
if isempty(regexp(mydir,'utility_scripts$','once'))
    warning('NO2_PROF_REPO assumes that it is within utility_scripts in the repo, this does not appear to be true (path = %s)',mydir)
end
repo_dir = regexprep(mydir, '[\\/]utility_scripts$', '');


end

