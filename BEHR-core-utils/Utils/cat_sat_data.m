function [ varargout ] = cat_sat_data( filepath, datafields, varargin )
%CAT_SAT_DATA( FILEPATH, DATAFIELDS ) Concatenates data from OMI .mat files
%   In some cases, one might wish to use satellite data from multiple days,
%   but we import OMI data and process BEHR data into daily files. This
%   function will load all the .mat files in the directory given by
%   FILEPATH and output a concatenated version of the data in the field or
%   fields given by DATAFIELDS, which should be a string or cell array of
%   strings. The output will be each requested data field as a separate
%   output, plus a cell array that gives the date, orbit, and pixel number
%   of each value in the first data field.
%
%   CAT_SAT_DATA( DATA, DATAFIELDS ) will concatenate all swaths in the
%   structure DATA for the fields specified in DATAFIELDS. If DATAFIELDS is
%   an empty array, all fields will be concatenated.
%
%   Parameter arguments are:
%
%       'prefix' - the prefix of the file names (the part before the date).
%       This function will use all files in the directory FILEPATH that
%       match the pattern <PREFIX>*.mat, replacing <PREFIX> with the given
%       prefix. This defaults to an empty string, i.e. by default all .mat
%       files will be matched.
%
%       'startdate' and 'enddate' - any valid representation of dates
%       (datenum or datestr). If only one is given, it will be treated as a
%       start or end point with the other unspecified (e.g. if you specify
%       startdate as 1-Jan-2015, this will operate on all files after
%       1-Jan-2015)
%
%       'newdim' - boolean, defaults to false. When true, each variable will
%       be concatenated along a new dimension (so a 2D variable will be
%       concatenated along the third dimension, a 3D one along the fourth).
%       When false, they will be concatenated in the along track dimension.
%
%       'vector' - boolean, defaults to false. When true, each swatch is
%       resized into a column vector before concatenation. Mutually
%       exclusive with 'newdim'.
%
%       'varname' - must be the string 'Data' or 'OMI'. 'Data' is default.
%       Indicates whether the structure concatenated should be the native
%       pixels ('Data') or the gridded pixels ('OMI').
%
%       'reject_args' - a cell array of arguments to be passed to
%       OMI_PIXEL_REJECT after the Data structure (so the second and later
%       arguments). If not given (or an empty array) then OMI_PIXEL_REJECT
%       is not called. If given, then pixels rejected by OMI_PIXEL_REJECT
%       will be either replaced with NaN or (if 'vector' is true) removed.
%       NOTE: in vector mode, all fields have rejected values removed.
%       However, in all other modes, only numerical fields have rejected
%       pixels set to NaN; e.g. logical fields will not be affected by
%       pixel filtering.
%
%       'struct_out' - if false (default) each concatenated field is
%       returned separately as its own output. If true, the fields are
%       returned in a scalar structure.
%
%       'DEBUG_LEVEL' - set to 0 to suppress debugging messages, defaults
%       to 1. Set to 'visual' to use the waitbar dialogue.
%
%   Josh Laughner <joshlaugh5@gmail.com> 10 Sept 2015

E=JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT PARSING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

if isstruct(filepath)
    Data = filepath;
    load_data = false;
else
    load_data = true;
    if ~ischar(filepath) || ~exist(filepath,'dir')
        E.badinput('filepath must be a string specifying a valid directory')
    end
end

if ischar(datafields)
    datafields = {datafields};
elseif isempty(datafields) && isstruct(filepath)
    datafields = fieldnames(Data);
elseif ~iscell(datafields) || any(~iscellcontents(datafields,'ischar'))
    E.badinput('datafields must be a string or cell array of strings')
end

p=inputParser;
p.addParameter('prefix','',@ischar);
p.addParameter('startdate',0);
p.addParameter('enddate','3000-01-01'); % pick a date far enough into the future that effectively all files will be before it
p.addParameter('newdim',false);
p.addParameter('vector',false);
p.addParameter('varname','Data');
p.addParameter('reject_args',{});
p.addParameter('struct_out', false);
p.addParameter('DEBUG_LEVEL',1,@(x) (ischar(x) || isnumeric(x) && isscalar(x)));

p.parse(varargin{:});
pout = p.Results;

prefix = pout.prefix;
startdate = pout.startdate;
enddate = pout.enddate;
newdim = pout.newdim;
vector_bool = pout.vector;
varname = pout.varname;
reject_args = pout.reject_args;
do_output_struct = pout.struct_out;
DEBUG_LEVEL = pout.DEBUG_LEVEL;

if ~ismember(varname,{'Data','OMI'})
    E.badinput('VARNAME must be either ''Data'' or ''OMI''');
end

if newdim && vector_bool
    E.badinput('NEWDIM and VECTOR are mutually exclusive')
end

wbbool = false;
if ischar(DEBUG_LEVEL)
    if strcmpi(DEBUG_LEVEL, 'visual')
        if isDisplay
            wbbool = true;
            DEBUG_LEVEL = 0;
        end
    else
        warning('Only the string ''visual'' for DEBUG_LEVEL will trigger the use of the waitbar.')
        DEBUG_LEVEL = 1;
    end
end

startdate = validate_date(startdate);
enddate = validate_date(enddate);
date_array = [startdate(:), enddate(:)];

if startdate > enddate
    E.badinput('startdate is later than enddate.')
end

if ~isscalar(newdim) || (~islogical(newdim) && ~isnumeric(newdim))
    E.badinput('The parameter newdim must be understood as a scalar logical.')
end

if ~isscalar(vector_bool) || (~islogical(vector_bool) && ~isnumeric(vector_bool))
    E.badinput('The parameter vector must be understood as a scalar logical.')
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% MAIN FUNCTION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
if load_data
    % Get all .mat files in the specified directory
    F = dir(fullfile(filepath, sprintf('%s*.mat',prefix)));
    
    if isempty(F)
        E.filenotfound('satellite .mat file');
    end
else
    F = 0;
end

% Prep output - we'll always include the date and swath as the last output
output_data = cell(1,numel(datafields)+1);

% Loop over all files (within the date limits given). Load the data
% variable, look for the datafields given, and add their data to the output
% which will be one long column vector.

if wbbool && load_data
    wb = waitbar(0,sprintf('Concatenating %s*.mat',strrep(prefix,'_','\_')));
elseif wbbool
    wb = waitbar(0,'Concatenating input structure');
end

for a=1:numel(F)
    if load_data
        [s,e] = regexp(F(a).name, '\d\d\d\d\d\d\d\d');
        filedate = datenum(F(a).name(s:e), 'yyyymmdd');
        in_date_ranges = filedate >= date_array(:,1) & filedate <= date_array(:,2);
        if ~any(in_date_ranges)
            continue
        end
        
        D = load(fullfile(filepath, F(a).name),varname);
        if ~isfield(D, varname)
            fprintf('%s does not contain the variable "%s", skipping\n', F(a).name, varname);
            continue
        else
            Data = D.(varname);
        end
    
        if DEBUG_LEVEL > 0
            fprintf('Loading file %s...\n',F(a).name);
        elseif wbbool
            waitbar(a/numel(F));
        end
    end
    
    for b=1:numel(datafields)
        if ~isfield(Data,datafields{b})
            if isstruct(F)
                E.callError('fieldnotfound','The field %s is not present in file %s',datafields{b},F(a).name);
            else
                E.callError('fieldnotfound','The field %s is not present in the input Data structure',datafields{b});
            end
        end
        
        for c=1:numel(Data)
            this_field = Data(c).(datafields{b});
            this_date_and_swath = format_date_swath_cell(Data(c), size(this_field));
            
            if ~isempty(reject_args)
                Data(c).Areaweight = ones(size(Data(c).Longitude));
                Data(c) = omi_pixel_reject(Data(c), reject_args{:});
                xx_keep = Data(c).Areaweight > 0;
                if vector_bool
                    this_field(~xx_keep) = [];
                    this_date_and_swath(~xx_keep) = [];
                else
                    % Logical fields cannot have values set to NaN. Since
                    % only numeric fields should really be filtered for low
                    % quality data, this shouldn't be a problem.
                    if isnumeric(this_field)
                        this_field(~xx_keep) = NaN;
                    end
                    % No need to set date and swath to NaNs, we just want
                    % to mark bad data as invalid.
                end
            end
            
            if newdim
                n = ndims(this_field);
                output_data{b} = cat(n+1, output_data{b}, this_field);
                if b == 1
                    output_data{end} = cat(n+1, output_data{end}, this_date_and_swath);
                end
            elseif vector_bool
                output_data{b} = cat(1, output_data{b}, this_field(:));
                if b == 1
                    output_data{end} = cat(1, output_data{end}, this_date_and_swath(:));
                end
            elseif ~newdim && ismatrix(Data(c).(datafields{b}))
                output_data{b} = cat(1, output_data{b}, this_field);
                if b == 1
                    output_data{end} = cat(1, output_data{end}, this_date_and_swath);
                end
            elseif ~newdim && ~ismatrix(Data(c).(datafields{b}))
                output_data{b} = cat(2, output_data{b}, this_field);
                if b == 1
                    output_data{end} = cat(2, output_data{end}, this_date_and_swath);
                end
            else
                E.notimplemented(sprintf('concat case: newdim = %d and ndims = %d',newdim,ndims(Data(c).(datafields{b}))));
            end
        end
    end
end

if wbbool
    close(wb);
end

if do_output_struct
    output_field_names = veccat(datafields, {'DateAndSwath'});
    % Have to put the date/orbit cell array into a parent cell array to
    % avoid creating a multielement structure 
    output_data{end} = {output_data{end}}; %#ok<CCAT1>
    varargout{1} = make_struct_from_field_values(output_field_names, output_data);
else
    varargout = output_data;
end

end

function ds_cell = format_date_swath_cell(data, sz)
E = JLLErrors;
if ~isscalar(data) || ~isstruct(data)
    E.badinput('DATA must be a scalar structure');
end

date = data.Date;
% Old files might have Swath as a matrix, and sometimes will incorrectly
% have a 0 in there, this should ensure we get the correct swath number as
% a scalar value
swath = unique(data.Swath(data.Swath > 0));

ds_cell = cell(sz);
for a=1:numel(ds_cell)
    ds_cell{a} = sprintf('%s_o%d_p%d',date,swath,a);
end
end