function [data,rLen,readOK,headerData]=readFORTRANRecord(fileID,dataStr,dTypeSize,headerNum,headerType,headerDTypeSize)
%READFORTRANRECORD Reads FORTRAN binary files
%   Reads unformatted sequential FORTRAN binary files already open in
%   MATLAB

% Check that the file ID is OK
if fileID < 0
    error('ReadFortranRecord:InvalidFile','Invalid file identifier');
end

if nargin<2
    dataStr = '*int32';
    dTypeSize=4;
end

headerData=[];

% Header data: sometimes there are strange data types at the start. If this
% is the case, use:
%   headerNum:      Number of foreign entries
%   headerType:     Data type of foreign entries
%   headerDTypeSize:Bytes per data type of foreign entries

% Get record length
rLen=fread(fileID,1,'*int32');

if strcmpi(dataStr,'seekpast')
    % Data not needed - ignore
    readOK=fseek(fileID,rLen+4,0); %#ok<NASGU>
    data=false;
    readOK=true;
else
    if exist('headerNum','var')
        adjustLen=headerNum*headerDTypeSize;
        headerData=fread(fileID,headerNum,headerType);
    else
        adjustLen=0;
    end

    % rLen is the number of bytes to read, take dTypeSize=number of bytes per
    % entry

    data=(fread(fileID,(rLen-adjustLen)/dTypeSize,dataStr));

    % Check that final integer matches starting integer
    ender=fread(fileID,1,'*int32');
    if ender~=rLen
        readOK=false;
    else
        readOK=true;
    end
end

end
