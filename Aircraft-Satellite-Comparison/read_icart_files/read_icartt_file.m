function [ Merge, header, data_table, field_bool ] = read_icartt_file( ict_file, field_bool )
%READ_ICARTT_FILE( ICT_FILE ) Read a single ICARTT file
%   [ Merge, header, data_table ] = read_icartt_file( ict_file ) Read in
%   the ICARTT formatted data file ICT_FILE and return the data in two
%   formats. "Merge" is a structure containing the data as fields,
%   "data_table" is a cell array of the CSV table in the ICARTT file, and
%   "header" contains the variable names for that table.

%%%%% INPUT PARSING %%%%%
E = JLLErrors;
DEBUG_LEVEL = 2;

if ~ischar(ict_file)
    E.badinput('ICT_FILE must be a string')
elseif ~exist(ict_file, 'file')
    E.badinput('Given ICT_FILE does not exist')
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%% MAIN FUNCTION %%%%%

% Create the data structure
Merge = struct;

Merge.metadata.file = ict_file;

fid = fopen(ict_file);
line1 = fgetl(fid); % Some merge files separate the two numbers on the first line with a comma, some don't. This gets the whole line...
N = sscanf(line1, '%d'); % ...then this gets the first number on that line.

if length(N) == 1; delim = ',';
else delim = '';
end

line_count = 2; % Keep track of what line we're on

% Skip lines 2 & 3
for l=2:3; fgetl(fid); line_count = line_count + 1; end

% Get the description of the file for the metadata
Merge.metadata.description = fgetl(fid); line_count = line_count + 1;
% Create the UTC unit field
Merge.Data.UTC.Fill = 'N/A';
Merge.Data.UTC.Unit = 'seconds after midnight in UTC time';

% Skip lines 5-6
for l=5:6; fgetl(fid); line_count = line_count + 1; end

% Get the data date.  Check that y contains at least 3 values, if it
% doesn't, the file probably uses space delimiters instead of commas.
line7 = fgetl(fid); line_count = line_count + 1;
y = sscanf(line7,'%d, %d, %d');
if length(y) < 3
    y = sscanf(line7,'%d %d %d');
end
curr_date = sprintf('%04d-%02d-%02d', y(1), y(2), y(3));
Merge.metadata.date = curr_date;

% Skip lines 8-9
for l=8:9; fgetl(fid); line_count = line_count + 1; end

% Get the number of fields we'll need to handle
numfields = str2double(fgetl(fid)); line_count = line_count + 1;
if ~exist('field_bool', 'var') || isempty(field_bool)
    field_bool = -1*ones(1,numfields+1); % the +1 accounts for the UTC field
elseif strcmpi(field_bool, 'existing')
    field_bool = 1*ones(1,numfields+1);
elseif strcmpi(field_bool, 'header')
    field_bool = zeros(1,numfields+1);
elseif ischar(field_bool)
    E.badinput('FIELD_BOOL (if given as a character) must be either ''existing'' or ''header''')
elseif isnumeric(field_bool) && numel(field_bool) ~= numfields+1
    E.badinput('FIELD_BOOL (if given as a numeric vector) must have a number of elements equal to the number of fields in this file (%d)', numfields+1);
elseif ~isnumeric(field_bool)
    E.badinput('FIELD_BOOL (if given) must be either one of the strings ''existing'' or ''header'' or a numeric vector')
end
% Skip one line
fgetl(fid); line_count = line_count + 1;

% Read in the fill values, creating the format spec appropriately based
% on the delimiter decided from the first line.
line12 = fgetl(fid);
if strcmpi(delim,',')
    fill_vals = sscanf(line12,repmat('%f, ',1,numfields));
else
    fill_vals = sscanf(line12,repmat('%f ',1,numfields));
end
line_count = line_count + 1;

% For each field, read the unit and fill value into the data structure
field_array = cell(1,numfields+1); field_array{1} = 'UTC';
for f = 1:numfields
    line = fgetl(fid); line_count = line_count + 1;
    comma_pos = strfind(line,',');
    if isempty(comma_pos) %Some files don't comma separate
        comma_pos = strfind(line,' ');
        if ~isempty(comma_pos)
            field = line(1:(comma_pos-1));
            unit = line((comma_pos+1):end);
        else
            % The TexAQS2000 data does not put units with the fields,
            % instead they just say "all in pptv except CO and CH4,
            % which are in ppbv". So those files are a pain.
            % Technically the UTC field will be TimeOpen_WAS, but since
            % I'm working with WAS rather than full merge file, won't
            % worry about that.
            field = line;
            unit = 'Unknown, see .ict file';
        end
    else
        field = line(1:(comma_pos-1));
        unit = strtrim(line((comma_pos+1):end));
    end
    
    % Sanitize the field name, characters that don't belong in variable
    % names will be removed. Also, prepend 'f_' to any field names that
    % start with a non-letter character (and so are invalid field names)
    field = regexprep(field,'\W','');
    if ~isempty(regexp(field(1),'[^a-zA-Z]')) %#ok<RGXP1> % regexp(field,'[^a-zA-Z]','ONCE') DOES NOT WORK: this line needs to test if the first character in 'field' is not a letter.
        field = ['f_',field]; 
    end 
    Merge.Data.(field).Unit = unit;
    Merge.Data.(field).Fill = fill_vals(f);
    
    field_array{f+1} = field;
end

% Read lines until we've hit the big table of values
while line_count < N(1)
    line = fgetl(fid);
    
    %Parse the line to see if it has any useful data. Some fields are given
    %specific names in the metadata; this is a hold over from earlier
    %versions of this code.
    if regcmp(line,'^\s*DATA_INFO')
        Merge.metadata.info = line(11:end);
    elseif regcmp(line,'^\s*ULOD_FLAG')
        value_string = regexp(line, '(?<=ULOD_FLAG:).*', 'match', 'once');
        Merge.metadata.upper_lod_flag = str2double(value_string);
    elseif regcmp(line, '^\s*ULOD_VALUE')
        value_string = regexp(line, '(?<=ULOD_VALUE:).*', 'match', 'once');
        Merge.metadata.upper_lod_value = value_string;
    elseif regcmp(line, '^\s*LLOD_FLAG')
        value_string = regexp(line, '(?<=LLOD_FLAG:).*', 'match', 'once');
        Merge.metadata.lower_lod_flag = str2double(value_string);
    elseif regcmp(line, '^\s*LLOD_VALUE')
        value_string = regexp(line, '(?<=LLOD_VALUE:).*', 'match', 'once');
        Merge.metadata.lower_lod_value = value_string;
    elseif regcmp(line, '^[\s\w\-]+:.*\w+')
        % If the line starts with some other field name and has at least
        % one character after the colon, add it as another field to the
        % metadata.
        
        % Split out the category and the value. If there's >1 colon, put
        % the extra ones back in
        split_line = strsplit(line, ':');
        field_name = regexprep(strtrim(split_line{1}),'\W','_');
        if ~regcmp(field_name, '^[a-zA-Z]')
            % Field names in structure must start with a letter
            field_name = sprintf('Field%s', field_name);
        end
        field_value = strjoin(split_line(2:end),':');
        Merge.metadata.(field_name) = strtrim(field_value);
    end
    
    line_count = line_count + 1;
end

% Read in the header row
header = cell(1,numfields+1);
manual_parse = 0;

for f = 1:numfields+1
    % Always read in the first time, then only continue scanning the
    % file if manual parsing is not required
    if ~manual_parse
        field = fscanf(fid,'%s,');
    end
    
    %Header rows with commas but not spaces present a problem for
    %fscanf.  If the field is too long, assume that it accidentally
    %read the whole line in instead of just one field, manually parse
    %the line and break out of the enclosing loop.
    %Only check this the first time; it'll either happen or it won't.
    if length(field)>20 && ~manual_parse && f == 1
        headerline = field;
        manual_parse = 1;
    end
    
    if manual_parse
        delim_pos = strfind(headerline,delim);
        if ~isempty(delim_pos)
            field = headerline(1:delim_pos(1)-1);
            headerline = headerline((delim_pos+1):end);
        else %if no more delimiters left, we're at the end of line (probably)
            field = deblank(headerline);
        end
    end
    
    field = regexprep(field,'\W','');
    if ~isempty(regexp(field(1),'[^a-zA-Z]')); field = ['f_',field]; end %#ok<RGXP1> % regexp(field,'[^a-zA-Z]','ONCE') DOES NOT WORK: this line needs to test if the first character in 'field' is not a letter.
    header{f} = field;
end

% Close the fid; we will use dlmread now and no longer need this file
% open
fclose(fid);

if DEBUG_LEVEL > 0; fprintf('   Loading the data table\n'); end
% Read in the full table
data_table = dlmread(ict_file,delim,N(1),0);

% Double check that the table size matches the header row
if size(data_table,2) ~= numel(header)
    error('read_merge:data_import','Table imported and header row do not have the same number of columns');
end

% Append the table data to the Merge structure
for r = 1:numfields+1
    new_field = header{r}; old_field = field_array{r};
    if field_bool(r) == 1
        field = old_field;
    elseif field_bool(r) == 0
        field = new_field;
    elseif field_bool(r) == -1
        if strcmp(new_field,old_field)
            field_bool(r) = 1;
            field = new_field;
        else
            question = sprintf('Header field %s does not match expected field %s. Which should be used?',new_field, old_field);
            choice = questdlg(question,'Field Mismatch','Header','Existing','Cancel run','Header');
            switch choice
                case 'Header'
                    field = new_field; field_bool(r) = 0;
                case 'Existing'
                    field = old_field; field_bool(r) = 1;
                case 'Cancel run'
                    error('read_merge:user','User cancelled run');
            end
        end
    end
    if DEBUG_LEVEL > 1; fprintf('\tField: %s', field); end
    data_column = data_table(:,r)';
    Merge.Data.(field).Values = data_column;
    if DEBUG_LEVEL > 1; fprintf(' Done.\n'); end
end

% Remove any fields without actual data and reorder them to match the order
% in the header
fns = fieldnames(Merge.Data);
for i_fn = 1:numel(fns)
    if ~isfield(Merge.Data.(fns{i_fn}), 'Values')
        Merge.Data = rmfield(Merge.Data, fns{i_fn});
    end
end

return

for i_field = 1:numel(header)
    reordered_data.(header{i_field}) = Merge.Data.(header{i_field});
end
Merge.Data = reordered_data;

end

