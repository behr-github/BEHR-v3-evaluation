function [ Data ] = read_geos_output(  )
%read_geos_output Uses the functions written by Sebastion Eastham to pull
%in variables from geos chem output (binary punch) files.
%  	The read BPCH function require a goodly number of inputs; this uses the
%  	MATLAB UI features to make it slightly more user friendly.  It will ask
%  	for the category (diagnostic) and tracer name - these should be the
%  	"Matlab sanitized" versions.
%
%   This function assumes that the tracerinfo.dat and diaginfo.dat files
%   are in the same folder as the output file of interest. If that is not
%   true, an error will occur.

% Initialize the error handling class
E = JLLErrors;

% Add the folder with the BPCH functions. In should be in the folder with
% this function.
addpath('~/Documents/MATLAB/GEOS_Chem_Utils/BPCH_Functions');

% Get the output file from the user, and use the directory to find the
% tracerinfo and diaginfo files.

if isDisplay
    getfile_title = 'Select the output file to read. tracerinfo and diaginfo should be in the same folder.';
    [output_filename, pathname] = uigetfile('*',getfile_title);
    if output_filename == 0;
        error(E.userCancel);
    end
else
    user_filepath = input('Enter the path to the bpch file', 's');
    if isempty(user_filepath)
        E.userCancel;
    elseif ~exist(user_filepath,'file')
        E.filenotfound('bpch')
    end
    [pathname, filename, fileext] = fileparts(user_filepath);
    output_filename = [filename, fileext];
end

if ~exist(fullfile(pathname,'tracerinfo.dat'),'file') || ~exist(fullfile(pathname,'diaginfo.dat'),'file')
    error(E.callError('dat_files_missing','The tracerinfo.dat and diaginfo.dat files must in the same folder as the output file.'));
end

inputFile = fullfile(pathname,output_filename);
tracerFile = fullfile(pathname,'tracerinfo.dat');
diagFile = fullfile(pathname,'diaginfo.dat');

Data = struct('dataBlock', [], 'dataUnit', [], 'fullName', [], 'fullCat', [], 'tVec', [], 'modelName', [], 'modelRes', [], 'dataScale', [], 'molMass', [], 'tEdge', []);
input_title = 'Enter the %s, or: \n\t"SAT" to return fields needed to calculate column density \n\tPROD to return fields of NO production \n\tCLOUDS to get cloud level and fraction \n\t Cancel to exit and return';

% Ask the user for the category and tracer names.  Loop until the user
% cancels or enters an empty string for category or tracer; this allows the
% user to enter multiple variables.  Note that inputdlg returns an empty
% cell array, i.e. {}, if the user cancels and a cell array with an empty
% string, i.e. {''}, if the user hits OK without entering anything.  We'll
% use this to check if the user wants to append the data useful for dealing
% with satellites, such as pressure, boxheight, etc.
first_time = true;
append_sat = false;
append_prod = false;
append_clouds = false;
while true
    if isDisplay
        user_input = inputdlg(sprintf(input_title,'category'));
    else
        temp_input = input(sprintf(input_title, 'category'), 's');
        user_input{1} = temp_input;
    end
    if isempty(user_input) % Canceling returns empty cell array
        break
    elseif isempty(user_input{1})
        % Hitting OK on an empty string before any tracers are appended is
        % probably a mistake - so let the user know and try again.
        if first_time
            if isDisplay
                box_msg = 'An empty entry does not return anything. If you want to exit with an empty structure, use "Cancel"';
                box_title = 'Empty category';
                uiwait(msgbox(box_msg,box_title,'error','modal'));
            else
                box_msg = '\nAn empty entry does not return anything. If you want to exit with an empty structure, press Ctrl+C\n';
                fprintf(box_msg);
            end
            continue
        else 
            break
        end
    elseif strcmpi(user_input{1},'sat')
        append_sat = true;
        break
    elseif strcmpi(user_input{1},'prod')
        append_prod = true;
        break
    elseif strcmpi(user_input{1},'clouds')
        append_clouds = true;
        break
    end
    
    category = user_input{1};
    
    if isDisplay
        user_input = inputdlg(sprintf(input_title,'tracer'));
    else
        temp_input = input(sprintf(input_title, 'tracer'), 's');
        user_input{1} = temp_input;
    end
    if isempty(user_input) % Canceling returns empty cell array
        break
    elseif isempty(user_input{1})
        % Hitting OK on an empty string before any tracers are appended is
        % probably a mistake - so let the user know and try again.
        if first_time
            if isDisplay
                box_msg = 'An empty entry does not return anything. If you want to exit with an empty structure, use "Cancel"';
                box_title = 'Empty tracer';
                uiwait(msgbox(box_msg,box_title,'error','modal'));
            else
                box_msg = '\nAn empty entry does not return anything. If you want to exit with an empty structure, press Ctrl+C\n';
                fprintf(box_msg);
            end
            continue
        else 
            break
        end
    elseif strcmpi(user_input{1},'sat')
        append_sat = true;
        break
    elseif strcmpi(user_input{1},'prod')
        append_prod = true;
        break
    elseif strcmpi(user_input{1},'clouds')
        append_clouds = true;
        break
    end
    
    tracer = user_input{1};
    
    % Change the message for successive variables.
    if first_time
        first_time = false;
        input_title = 'Enter another %s, or: \n\t"SAT" to append fields needed to calculate column density \n\tPROD to append fields of NO production \n\tCLOUDS to get cloud level and fraction \n\t Press cancel or enter an empty string to exit and return.';
    end
    
    % Call the function to read the BPCH file
    try
        Data = read_tracer(category, tracer, inputFile, tracerFile, diagFile, Data);
    catch err
        % Catch errors resulting from the user entering a wrong category or
        % tracer.  This will present a message box that waits until the
        % user presses OK before continuing.
        if strcmp(err.identifier,'BPCHRead:UnknownTracer')
            tr_or_cat = 'tracer';
            couldnt_find = tracer;
        elseif strcmp(err.identifier,'BPCHRead:UnknownCategory')
            tr_or_cat = 'category';
            couldnt_find = category;
        else
            rethrow(err);
        end
        box_msg = sprintf('The %s %s was invalid. This category & tracer will not be read.',tr_or_cat,couldnt_find);
        box_title = sprintf('Could not find %s',tr_or_cat); 
        uiwait(msgbox(box_msg,box_title,'error','modal'));
    end
end

% Get the NAIR, PSURF, and BXHEIGHT variables, if the user hit "OK" on a
% blank input
if append_sat
    categories = {'BXHGHT','PEDGE','BXHGHT','TR_PAUSE'};
    tracers = {'NAIR','PSURF','BXHEIGHT','TP_LEVEL'};
    for a=1:numel(categories);
        Data = read_tracer(categories{a}, tracers{a}, inputFile, tracerFile, diagFile, Data);
    end
end

% Get all 8 categories for production of NO
if append_prod
    categories = {'NO_AC','NO_AN','NO_BIOB','NO_BIOF','NO_FERT','NO_LI','NO_SOIL','NO_STRT'};
    tracers = {'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO', 'NO'};
    for a=1:numel(categories)
        Data = read_tracer(categories{a}, tracers{a}, inputFile, tracerFile, diagFile, Data);
    end
end

% Get cloud top pressure and cloud fraction
if append_clouds
    categories = {'DAO_FLDS','DAO_FLDS'};
    tracers = {'CLDFRC','CLDTOP'};
    for a=1:numel(categories)
        Data = read_tracer(categories{a}, tracers{a}, inputFile, tracerFile, diagFile, Data);
    end
end

end

%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% SUBFUNCTIONS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function Data = read_tracer(category, tracer, inputFile, tracerFile, diagFile, Data)
    fprintf('Now retrieving %s/%s\n',category,tracer);
    [ dataBlock, dataUnit, fullName, fullCat, tVec, modelName, modelRes, dataScale, molMass, tEdge ] = readBPCHSingle(inputFile,category,tracer,tracerFile,diagFile);
    D = numel(Data);
    if ~isempty(Data(D).dataBlock)
        % Initializing Data as an empty structure puts empty matrices in
        % each field, but still makes it be a 1x1 structure. This check
        % figures out if its 1x1 because its freshly initialized, or
        % because only one tracer was read in.
        D=D+1;
    end
    Data(D).dataBlock = dataBlock;
    Data(D).dataUnit = dataUnit;
    Data(D).fullName = fullName;
    Data(D).fullCat = fullCat;
    Data(D).tVec = tVec;
    Data(D).modelName = modelName;
    Data(D).modelRes = modelRes;
    Data(D).dataScale = dataScale;
    Data(D).molMass = molMass;
    Data(D).tEdge = tEdge;
end
