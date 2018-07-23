classdef GitChecker < handle
    %GitChecker Use to verify that a dependency's Git repo is in the proper state
    %   Particularly in Matlab, where you can set a very wide search path
    %   to look for .m files, you may have code in one or more Git
    %   repositories that depends on code in another. Because Matlab will
    %   look through its search path and find the first file that matches a
    %   function call, there's no (easy) way to have two versions of
    %   function foo() in repository D that are needed for functions bar()
    %   and baz() in repositories A and B. If foo() has two different
    %   versions in two branches of D, you run the risk of using the wrong
    %   version of foo() for bar() or baz().
    %
    %   This class offers a way to prevent that. By adding an instance of
    %   this class to the beginning of a function call, you can verify that
    %   a second repository is in the state you want.
    %
    %   Usage: create an instance of this class, e.g. G = GitChecker(), and
    %   use the methods addAllowedBranches(), addReqCommits(), and
    %   addCommitRange() to specify rules for one or more repositories.
    %   (Use the classhelp function in this repository to get more
    %   information on those methods.)
    %
    %   There are also two properties which can be set to true to control
    %   the behavior of this instance. Strict will, if true, cause an error
    %   to automatically be thrown if any of the repository state
    %   conditions are not met (this is the default).  You can set this to
    %   false instead so that the check methods return a boolean for
    %   success or failure, which you must then handle yourself. Verbose
    %   will cause more information to be printed to the screen if set to
    %   true (default is false).
    
    properties
        Strict = true;
        Verbose = false;
    end
    properties(SetAccess = protected)
        Dirs = struct('directory','','allowed_branches',{{}},'req_commits',{{}},'commit_ranges',{{}})
    end
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% FUNCTIONS FOR ADDING RULES %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function addAllowedBranches(obj, gitdir, varargin)
            % G.addAllowedBranches(gitdir, branch1, branch2, ...)
            %   Add a rule saying that the repo GITDIR must be on one of
            %   the specified branches. The branches can either all be
            %   specified is a single call or by successive calls. GITDIR
            %   must be a path to any folder in the desired repo (the top
            %   level repo directory is fine)
            ind = obj.findOrAddDir(gitdir);
            for a=1:numel(varargin)
                branchname = varargin{a};
                if obj.isValidBranch(gitdir,branchname)
                    obj.Dirs(ind).allowed_branches{end+1} = branchname;
                else
                    error('gitchecker:branch_name','%s is not a valid branch for the repo in %s', branchname, gitdir);
                end
            end
        end
        
        function addReqCommits(obj, gitdir, varargin)
            % G.addReqCommitsByHash(gitdir, hash1, hash2, ...)
            %   Add a rule requiring that all the given commits are
            %   ancestors of the current HEAD. It is best if you identify
            %   the commits by their hash or tag rather than a branch name
            %   or a relative commit (i.e. HEAD~4) because those can change
            %   over time. The commits can either all be specified is a
            %   single call or by successive calls. GITDIR must be a path
            %   to any folder in the desired repo (the top level repo
            %   directory is fine)
            ind = obj.findOrAddDir(gitdir);
            for a = 1:numel(varargin)
                hash = varargin{a};
                if obj.isValidHash(gitdir, hash)
                    obj.Dirs(ind).req_commits{end+1} = hash;
                else
                    error('gitchecker:hash_id','%s is not a valid commit ID for the repo in %s', hash, gitdir);
                end
            end
        end
        
        function addCommitRange(obj, gitdir, first_hash, last_hash)
            % G.addCommitRangeByHash(gitdir, first_hash, last_hash)
            %   Add a rule requiring that the current HEAD is between
            %   FIRST_HASH and LAST_HASH (inclusive) in the repo history,
            %   or more specifically, that FIRST_HASH is an ancestor of
            %   HEAD and HEAD is an ancestor of LAST_HASH. Unlike
            %   addAllowedBranches and addReqCommits, this can only take
            %   one range at a time. If multiple ranges are specified, the
            %   HEAD can be in any of them. GITDIR must be a path to any
            %   folder in the desired repo (the top level repo directory is
            %   fine)
            ind = obj.findOrAddDir(gitdir);
            
            if ~ischar(first_hash) || ~ischar(last_hash)
                error('gitchecker:bad_input','Both FIRST_HASH and LAST_HASH must be strings');
            elseif ~obj.isValidHash(gitdir, first_hash)
                error('gitchecker:bad_input','FIRST_HASH (%s) does not appear to be a valid hash for repo in %s', first_hash, gitdir)
            elseif ~obj.isValidHash(gitdir, last_hash)
                error('gitchecker:bad_input','LAST_HASH (%s) does not appear to be a valid hash for repo in %s', last_hash, gitdir)
            end
            
            obj.Dirs(ind).commit_ranges{end+1} = {first_hash, last_hash};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%% CHECKING FUNCTIONS %%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function stat = checkState(obj, gitdir)
            if ~exist('gitdir', 'var')
                gitdir = '';
            end
            
            stat = obj.checkAllowedBranches(gitdir);
            stat = stat && obj.checkReqCommits(gitdir);
            stat = stat && obj.checkCommitRanges(gitdir);
            
        end
        
        function stat = checkAllowedBranches(obj, gitdir)
            if ~exist('gitdir', 'var')
                gitdir = '';
            end
            inds = obj.checkFxnIndicies(gitdir);
            stats = true(size(inds));
            currdir = pwd;
            for a=1:numel(inds)
                i = inds(a);cd(obj.Dirs(i).directory);
                try % encase in try-catch to always be sure to return to old directory in case of error
                    for b=1:numel(obj.Dirs(i).allowed_branches)
                        sys_cmd = sprintf('git rev-parse %s', obj.Dirs(i).allowed_branches{b});
                        [cmd_stat, branch_hash] = system(sys_cmd);
                        if cmd_stat
                            error('gitchecker:check_branches', 'Error during system call to obtain commit hash for branch %s', obj.Dirs(i).allowed_branches{b});
                        end
                        [cmd_stat, head_hash] = system('git rev-parse HEAD');
                        if cmd_stat
                            error('gitchecker:check_branches', 'Error during system call to obtain commit hash for HEAD');
                        end
                        stats(a) = strcmp(branch_hash, head_hash);
                        if stats(a)
                            % If on any of the allowed branches, this check
                            % passes. No need to do the other checks.
                            break
                        end
                    end
                catch err
                    cd(currdir);
                    rethrow(err);
                end
                cd(currdir);
                
                % If we get here, the test should have failed, but let's be
                % sure
                if ~stats(a)
                    msg = sprintf('Current HEAD for %s not on any of the allowed branches: %s', obj.Dirs(i).directory, strjoin(obj.Dirs(i).allowed_branches, ', '));
                    if obj.Strict
                        error('gitchecker:checkAllowedBranches', msg) %#ok<SPERR>
                    elseif obj.Verbose
                        fprintf('%s\n', msg);
                    end
                end
            end
            
            stat = all(stats);   
        end
        
        function stat = checkReqCommits(obj, gitdir)
            if ~exist('gitdir', 'var')
                gitdir = '';
            end
            inds = obj.checkFxnIndicies(gitdir);
            stats = true(size(inds));
            currdir = pwd;
            for a=1:numel(inds)
                i = inds(a);
                cd(obj.Dirs(i).directory);
                try % encase in try-catch to always be sure to return to old directory in case of error
                    for b=1:numel(obj.Dirs(i).req_commits)
                        sys_cmd = sprintf('git merge-base --is-ancestor %s HEAD', obj.Dirs(i).req_commits{b});
                        stats(a) = ~system(sys_cmd);
                        if ~stats(a)
                            % If any check fails, all fail, so no need to
                            % continue checking
                            continue
                        end
                    end
                catch err
                    cd(currdir);
                    rethrow(err);
                end
                cd(currdir);
                
                if ~stats(a)
                    msg = sprintf('Required commits (%s) not ancestors of current HEAD for %s\n', strjoin(obj.Dirs(i).req_commits, ', '), obj.Dirs(i).directory);
                    if obj.Strict
                        error('gitchecker:checkReqCommits', msg); %#ok<SPERR>
                    elseif obj.Verbose
                        fprintf('%s\n', msg);
                    end
                end
            end
            stat = all(stats);
        end
        
        function stat = checkCommitRanges(obj, gitdir)
            if ~exist('gitdir', 'var')
                gitdir = '';
            end
            inds = obj.checkFxnIndicies(gitdir);
            stats = true(size(inds));
            currdir = pwd;
            for a=1:numel(inds)
                i = inds(a);
                cd(obj.Dirs(i).directory);
                try % encase in try-catch to always be sure to return to old directory in case of error
                    for b=1:numel(obj.Dirs(i).commit_ranges)
                        sys_cmd = sprintf('git merge-base --is-ancestor %s HEAD', obj.Dirs(i).commit_ranges{b}{1});
                        stat1 = ~system(sys_cmd);
                        sys_cmd = sprintf('git merge-base --is-ancestor HEAD %s', obj.Dirs(i).commit_ranges{b}{2});
                        stat2 = ~system(sys_cmd);
                        stats(a) = stat1 && stat2;
                        if stats(a)
                            % If we are in at least one of the ranges, this
                            % test passes.
                            continue
                        end
                    end
                catch err
                    cd(currdir);
                    rethrow(err);
                end
                cd(currdir);
                
                if ~stats(a)
                    msg = sprintf('HEAD for %s not in any of the allowed commit ranges', obj.Dirs(i).directory);
                    if obj.Strict
                        error('gitchecker:checkCommitRanges', msg); %#ok<SPERR>
                    elseif obj.Verbose
                        fprintf('%s\n', msg);
                    end
                end
            end
            stat = all(stats);
        end
    end
    
    methods(Access=protected)
        function ind = findOrAddDir(obj, newdir, noadd)
            % Give the right index to add the branch or commit to in the
            % Dirs property. If the new directory isn't listed in Dirs, add
            % it. Also, always refer to the root directory of the
            % repository.
            if ~ischar(newdir)
                error('gitchecker:bad_input','Given GITDIR must be a string');
            elseif ~exist(newdir, 'dir')
                error('gitchecker:bad_input','Given GITDIR is not an extant directory');
            end
            
            if ~exist('noadd', 'var')
                noadd = false;
            elseif ~isscalar(noadd) || ~islogical(noadd)
                error('gitchecker:bad_input','NOADD must be a scalar logical value, or be omitted')
            end
            
            gitdir = obj.getGitRoot(newdir);
            extant_dirs = {obj.Dirs.directory};
            xx = strcmp(extant_dirs, gitdir);
            if sum(xx) > 0
                ind = find(xx);
                return
            end
            
            % Should really only be invoked the first time when we have an
            % initialized structure with an empty string for the directory.
            xx = strcmp(extant_dirs,'');
            if sum(xx) > 0
                ind = find(xx,1);
                obj.Dirs(ind).directory = gitdir;
                return
            end
            
            % Last case - doesn't exist and no empty directories present
            if noadd
                error('gitchecker:bad_git_dir', 'Git directory %s has not been specified in any existing rules', gitdir);
            else
                obj.Dirs(end+1).directory = gitdir;
                ind = numel(obj.Dirs);
            end
        end
        
        function inds = checkFxnIndicies(obj, gitdir)
            if isempty(gitdir)
                inds = 1:numel(obj.Dirs);
            elseif ischar(gitdir)
                inds = obj.findOrAddDir(gitdir, true);
            elseif isnumeric
                if any(gitdir < 1 | gitdir > numel(obj.Dirs))
                    error('gitchecker:bad_git_dir', 'If given as an array of indices, GITDIR must have all values between 1 and numel(obj.Dirs), currently %d', numel(obj.Dirs));
                else
                    inds = unique(gitdir);
                end
            end
        end
    end
    
    methods(Static)
        function p = getAbsolutePath(p)
            currdir = cd(p);
            p = pwd;
            cd(currdir);
        end
        
        function gitroot = getGitRoot(p)
            currdir = cd(p);
            syscmd = 'git rev-parse --show-toplevel';
            [stat,gitroot] = system(syscmd);
            cd(currdir);
            if stat == 128
                error('gitchecker:not_repo','%s is not in a git repository',p);
            elseif stat ~= 0
                error('gitchecker:sys_call','Error with system call "%s": %s', syscmd, gitroot)
            end
            gitroot = strtrim(gitroot);
        end
        
        function bool = isValidBranch(gitdir, branchname)
            syscmd = sprintf('git show-ref --verify --quiet refs/heads/%s',branchname);
            currdir = cd(gitdir);
            stat = system(syscmd);
            cd(currdir)
            bool = ~stat;
        end
        
        function bool = isValidHash(gitdir, hash)
            syscmd = sprintf('git rev-parse --quiet --verify %s',hash);
            currdir = cd(gitdir);
            [stat, ~] = system(syscmd);
            cd(currdir)
            bool = ~stat;
        end
    end
    
end

