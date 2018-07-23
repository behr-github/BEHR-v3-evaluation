function [ githead ] = git_head_hash( git_dir )
%GIT_HEAD_HASH Return the hash of the HEAD of a Git repository
%   GITHEAD = GIT_HEAD_HASH( GIT_DIR ) Returns the full hexadecimal SHA
%   commit hash of the commit pointed to by HEAD in the repository GIT_DIR.
%   Requires that a system call to "git" is possible on this machine.

currdir = cd(git_dir);
try
    [gitstat, githead] = system('git rev-parse HEAD');
catch err
    cd(currdir);
    rethrow(err);
end
cd(currdir);

if gitstat ~= 0
    githead = 'Unknown';
end

githead = strtrim(githead);

end

