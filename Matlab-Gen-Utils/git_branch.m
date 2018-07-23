function [ gbranch ] = git_branch( git_dir )
%GIT_BRANCH Get the current Git branch for a given directory
%   GBRANCH = GIT_BRANCH( GIT_DIR ) Returns the current branch for the
%   repository in GIT_DIR as a string. If the git command fails for any
%   reason, GBRANCH will be "Unknown".

olddir = cd(git_dir);
try
    [git_stat, gbranch] = system('git rev-parse --abbrev-ref HEAD');
catch err
    cd(olddir);
    rethrow(err);
end
cd(olddir);

if git_stat ~= 0
    gbranch = 'Unknown';
end

gbranch = strtrim(gbranch);

end

