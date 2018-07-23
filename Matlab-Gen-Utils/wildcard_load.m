function [ S ] = wildcard_load( filedir, filepattern, date_in, varargin )
%wildcard_load - Allows the user to load a file using wildcards
%   The built-in Matlab load command requires a totally specified file name
%   and path; but sometimes a file (especially a satellite file) can be
%   exactly specified with only a part of the file name.  Since wild cards
%   can only be used to get files through the dir command, this function
%   uses that to find the file desired.
%
%   The loaded variables are returned as a structure (the same structure
%   that is returned from the built in load command).  This seems to be the
%   best compromise between making clear what variables are returned and
%   minimizing the variable reassignment the user must do outside this
%   function.
%
%   This will automatically fill in dates in a file pattern if the third
%   argument is a date (number or valid string - see next paragraph for
%   necessary formatting of file pattern).  If an empty matrix (i.e. [] )
%   is passed as the third argument, this function will not try to fill in
%   a date.
%
%   User inputs: the file directory, the file pattern, the date for
%   the file (as a datenum), and (optionally) the variables to load.  The
%   file pattern must contain three '%s' format specs used by sprintf; the
%   first one will be the 4 digit year, the second the 2 digit month and
%   the third the 2 digit day.  If the year, month, day do not occur in
%   this order, use '%1$s','%2$s', and '%3$s' respectively.

    E = JLLErrors;
    
    
    if ~isempty(date_in)
        curr_date = datestr(date_in,29);
        year = curr_date(1:4);
        month = curr_date(6:7);
        day = curr_date(9:10);
        filename = sprintf(filepattern,year,month,day);
    else
        filename = filepattern;
    end
    
    file_list = dir(fullfile(filedir,filename));
    if numel(file_list)==1
        S = load(fullfile(filedir, file_list(1).name),varargin{:});
    elseif isempty(file_list)
        S = struct('a',{});
    else
        error(E.toomanyfiles);
    end


end

