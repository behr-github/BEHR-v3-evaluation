function [ builtin_bool ] = isbuiltin( class_in )
%ISBUILTIN Returns true if a function/class is part of MATLAB or a toolbox
%   Checks the path of a function or class for the usual location of
%   built in MATLAB functions or toolboxes. Returns true if it finds it,
%   indicating that the function shipped with MATLAB. One caveat: if you're
%   putting your own functions in the directory returned by matlabroot,
%   this will get confused.
%
%   Josh Laughner <joshlaugh5@gmail.com> 17 June 2015

C = which(class_in);
r = regexp(C,matlabroot,'ONCE');
if isempty(r)
    builtin_bool = false;
else
    builtin_bool = true;
end

end

