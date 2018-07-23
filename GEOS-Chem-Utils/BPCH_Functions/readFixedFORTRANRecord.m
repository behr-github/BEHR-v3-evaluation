function [readOK,varargout]=readFixedFORTRANRecord(fileID,varargin)
%READFIXEDFORTRANRECORD Reads FORTRAN binary files
%   Reads unformatted sequential FORTRAN binary files already open in
%   MATLAB. Expects arguments in pairs following the file ID, in groups
%   as follows:
%       dType:  String. Datatype, e.g. '*char'
%       dNum:   Integer. Number of variables in section.
%   Data will be returned, in order, as varargout.

% Get record length
rLen=fread(fileID,1,'*int32');

% Determine number of records to generate
numRecords = length(varargin)/2;
varargout = cell(numRecords,1);

for iRec = 1:numRecords
    varargout(iRec) = {fread(fileID,varargin{2*iRec},char(varargin{(2*iRec)-1}))};
end

% Check that final integer matches starting integer
ender=fread(fileID,1,'*int32');
if ender~=rLen
    readOK=false;
else
    readOK=true;
end

end
