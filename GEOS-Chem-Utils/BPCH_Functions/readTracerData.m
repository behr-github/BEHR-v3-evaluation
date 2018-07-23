function [tID,tName,tWeight,tCarbon,tNum,tScale,tUnit]=readTracerData(tracerInfo)

% Use persistent hash table to save time
persistent tracerMap
dataStored = false;
makeHashTable = false;
if isempty(tracerMap)
    makeHashTable = true;
elseif tracerMap.isKey(tracerInfo)
    dataStored = true;
end

if ~dataStored
    % Read tracer info - expect formatted tracerinfo file
    tinfoID = fopen(tracerInfo,'r','ieee-be');
    % Scan file to determine array sizes
    % Would be slightly faster to simply declare large arrays first and then
    % trim, but this is more memory efficient
    tracerBlock = false;
    tCount = 0;
    while ~feof(tinfoID)
        if ~tracerBlock
            startPt = ftell(tinfoID);
        end
        currLine = fgetl(tinfoID);
        if ~strcmpi('#',currLine(1));
            tCount = tCount + 1;
            tracerBlock = true;
        end
    end

    % Declare relevant arrays
    tID     = cell(tCount,1);
    tName   = cell(tCount,1);
    tWeight = zeros(tCount,1);
    tCarbon = zeros(tCount,1);
    tNum    = zeros(tCount,1);
    tScale  = zeros(tCount,1);
    tUnit   = cell(tCount,1);

    % Rewind to the start of the tracer block
    fseek(tinfoID,startPt,-1);

    % Read in tracers
    iTracer = 0;
    while ~feof(tinfoID)
        currLine = fgetl(tinfoID);
        if ~strcmpi('#',currLine(1));
            % Valid line
            iTracer = iTracer + 1;
            tID(iTracer)    = {strtrim(currLine(1:8))};
            tName(iTracer)  = {strtrim(currLine(10:40))};
            tWeight(iTracer)= str2double(currLine(41:49));
            tCarbon(iTracer)= str2double(currLine(50:52));
            tNum(iTracer)   = str2double(currLine(53:61));
            tScale(iTracer) = str2double(currLine(62:71));
            tUnit(iTracer)  = {strtrim(currLine(72:end))};
        end
    end
    fclose(tinfoID);
    % Store data in hash table
    if makeHashTable
        % Establish hash table
        tracerMap = containers.Map({tracerInfo},{{tID,tName,tWeight,tCarbon,tNum,tScale,tUnit}});
    else
        % Just add key
        tracerMap(tracerInfo) = {tID,tName,tWeight,tCarbon,tNum,tScale,tUnit};
    end
else
    % Retrieve from hash table
    tracerData = tracerMap(tracerInfo);
    tID = tracerData{1};
    tName = tracerData{2};
    tWeight = tracerData{3};
    tCarbon = tracerData{4};
    tNum = tracerData{5};
    tScale = tracerData{6};
    tUnit = tracerData{7};
end
end