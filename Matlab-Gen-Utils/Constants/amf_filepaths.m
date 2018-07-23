function [ fileDamf, fileTmp ] = amf_filepaths( )
%amf_filepaths Returns the file paths for the dAmf file (scattering weights
%table) and the nmcTmpYr file (temperature profiles for correction factor
%alpha)
%   Update as needed if file paths change.  This function should always be
%   referenced for these paths to make code as portable as possible.

global onCluster
if isempty(onCluster)
    onCluster = 0;
end

if onCluster
    amf_tools_path = '/global/home/users/laughner/MATLAB/BEHR/AMF_tools';
else
    homedir = getenv('HOME');
    amf_tools_path = sprintf('%s/Documents/MATLAB/BEHR/AMF_tools',homedir);
end

fileTmp = fullfile(amf_tools_path,'nmcTmpYr.txt');
fileDamf = fullfile(amf_tools_path,'damf.txt');

end

