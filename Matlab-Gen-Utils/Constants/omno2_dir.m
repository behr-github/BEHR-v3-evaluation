function [ omno2_path ] = omno2_dir( age )
%omno2_dir Returns the path to the OMNO2 files. 
%   Returns a string with the path to the OMNO2 files.  Defaults to the old
%   file server, pass the string 'new' to use the new one.

if nargin > 0 && strcmpi(age,'old');
    omno2_path = '/Volumes/share/GROUP/SAT/OMI/OMNO2_32';
else
    omno2_path = '/Volumes/share-sat/SAT/OMI/OMNO2';
end



end

