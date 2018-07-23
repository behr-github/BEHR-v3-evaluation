function [dOffset,dName,dFull]=readDiagInfo(diagInfo)
% Read diagnostic info - expect formatted diaginfo file
dinfoID = fopen(diagInfo,'r','ieee-be');
% Scan file to determine total number of diagnostics
diagBlock = false;
dCount = 0;
while ~feof(dinfoID)
    if ~diagBlock
        startPt = ftell(dinfoID);
    end
    currLine = fgetl(dinfoID);
    if ~strcmpi('#',currLine(1));
        dCount = dCount + 1;
        diagBlock = true;
    end
end

% Declare relevant arrays
dOffset = zeros(dCount,1);
dName   = cell(dCount,1);
dFull   = dName;

% Rewind to the start of the tracer block
fseek(dinfoID,startPt,-1);

% Read in tracers
iDiag = 0;
while ~feof(dinfoID)
    currLine = fgetl(dinfoID);
    if ~strcmpi('#',currLine(1));
        % Valid line
        iDiag = iDiag + 1;
        dOffset(iDiag)  = str2double(currLine(1:8));
        if length(currLine) > 49
            dName(iDiag)    = {strtrim(currLine(10:49))};
            dFull(iDiag)    = {strtrim(currLine(50:end))};
        else
            dName(iDiag)    = {strtrim(currLine(10:end))};
            dFull(iDiag)    = {strtrim(currLine(10:end))};
        end
    end
end
fclose(dinfoID);
end