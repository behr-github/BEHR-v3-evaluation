function [ behr_path ] = behr_dir(  )
%BEHR_DIR Returns the mounted path for BEHR files
%   Detailed explanation goes here

if ismac
    remote_path = '/Volumes';
else
    error('BEHR_DIR remote path not defined for your OS')
end

behr_path = fullfile(remote_path,'share-sat','SAT','BEHR','BEHR_Files_2014');

end

