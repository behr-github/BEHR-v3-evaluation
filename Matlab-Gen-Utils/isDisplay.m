function [ display_bool ] = isDisplay(  )
%isDisplay Try to determine if matlab started with the -nodisplay parameter
%   Matlab can be started from the command line with the -nodisplay option.
%   This disables certain features, such as dialogue boxes and the wait
%   bar. Since there is no good way to actually test if Matlab has a
%   display open ( get(0,'screensize') does not seem to reliably return 0's
%   for the last two arguments when -nodisplay is set, as some forum posts
%   suggest), I came up with a workaround.  
%
%   In my .bashrc file, I define an alias "startmatlab" that starts Matlab
%   in a terminal with -nodisplay and -nosplash.  In that alias, I also
%   define an environmental variable, MATLAB_DISPLAY, which is set to 0.
%   This function checks for that variable and returns a 1 or 0 to say if
%   the display is set.
%
%   Note that this takes advantage of how *nix systems (Mac and Linux)
%   work.  There may be no way to duplicate this on a PC.
%
%   Josh Laughner <joshlaugh5@gmail.com> 27 Feb 2015

home = getenv('HOME');


warnme = true; % Set this to false to disable the following warning
if ~strcmp(home,'/Users/Josh') && warnme
    warning('If you do not define a MATLAB_DISPLAY environmental variable in a shell when launching MATLAB with the -nodisplay option, this function will always return true. See documentation.\nThis warning can be disabled by setting "warnme = false" in the function');
end


display_bool = str2double(getenv('MATLAB_DISPLAY'));
if isnan(display_bool);
    display_bool = 1;
end

end

