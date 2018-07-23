function [ varargout ] = make_latex_table( T, varargin )
%MAKE_LATEX_TABLE Create a Latex table environment
%   TABLE_STR = MAKE_LATEX_TABLE( T ) Returns a string representation of
%   T that includes \begin{table} and \begin{tabular} (and corresponding
%   \end commands) environments. T may be an array, a table, or a string.
%   If given as a string, it must contain newlines processed by
%   sprintf or fprintf to transform the "\n" marker to the actual newline
%   character. If given as a table or array, it is passed to mat2latex to
%   generate the required string. If given as a table, the table
%   VariableNames are used as the column headers and the RowNames as the
%   row headers.
%
%   This function has a large number of parameter arguments that alter its
%   function:
%
%       'file' - allows the string to be directly output to a file.
%       Normally it just writes the string alone to that file, however
%       'insert' changes this behavior (see below). By default this will
%       ask before overwriting an existing file.
%
%       'insert' - instead of writing to a new file, this will instead try
%       to insert the table string into an existing file. It looks for
%       lines '%BEGIN <mark>' and '%END <mark>' and puts the table in
%       between those lines. <mark> is given by the 'marker' parameter (see
%       below). The %BEGIN or %END must be at the beginning of the line and
%       followed by at least one space before the <mark>. Everything
%       currently between the %BEGIN and %END lines is lost. This will back
%       up the existing file by copying it to <file>.bckp - but be wary of
%       running this function multiple times, as each time it will
%       overwrite the previous backup.
%
%       'marker' - the identifying string in the Latex file that it places
%       the table between in 'insert' mode. That is, if the value for this
%       is 'TABLE1', the table is placed between the lines "%BEGIN TABLE1"
%       and "%END TABLE1".
%
%       'overwrite' - boolean, if true, immediately overwrites the output
%       file without asking IF not operating in "insert" mode.
%
%       'caption' - the string to be placed in the \caption{} command for
%       the table environment. If not given, or if an empty string, the
%       \caption{} command will not be inserted at all.
%
%       'label' - the string to be place in the \label{} command for the
%       table environment. If not given, or if an empty string, the
%       \caption{} command will not be inserted at all.
%
%       'rownames' - cell array of strings that will be used as the row
%       names; overrides those given in T if it is a table. rownames must
%       by R-by-N, where R is the number of rows in the table and N is >=
%       1. By default, vertically adjacent identical strings are combined
%       into \multirow commands (see next parameter).
%
%       'multirow' - boolean, if true (default), then identical adjacent
%       row names along the first dimension will be combined into a
%       \multirow command. E.g. if 'rownames' was:
%
%           {   'alpha', 'ex'    ;
%               'alpha', 'ex'    ;
%               'alpha', 'why'   ;
%               'beta',  'zee'   ;
%               'alpha', 'zee'   }
%
%       then the first three "alpha"s and the two "ex"s and "zee"s would be
%       combined into \multirow{3}{*}{alpha}, \multirow{2}{*}{ex} and
%       \multirow{2}{*}{zee} respectively. However, the last "alpha" would
%       stay separate. If you give 'false' to this parameter, then no
%       entries will be combined into \multirow.
%
%       'colnames' - cell array of strings that will be used as the column
%       names; overrides the VariableNames given in T if it is a table.
%       colnames must be N-by-C where N is >= 1 and C is either the number
%       of columns in the table or the number of columns in the table plus
%       the number of columns in rownames. By default, horizontally
%       adjacent identical strings are combined into \multicolumn commands
%       (see next parameter).
%
%       'multicol' - boolean, if true (default), combines adjacent
%       identical column names into \multicolumn commands. This checks
%       along the second dimension of colnames, so e.g.:
%
%           { 'alpha', 'alpha', 'alpha' ;
%             'ex',    'why',   'zee'   }
%
%       would combine the three "alpha" into a \multicolumn{3}{c}{alpha}
%       command. Currently, there is no way to change the "c" format.
%       Setting this parameter to false means that nothing would be
%       combined.
%
%       'm2l' - a cell array of options that mat2latex understands. These
%       will be passed to mat2latex if mat2latex is called to convert T to
%       a string.
%
%       'lines' - a 1-by-3 cell array of strings that have Latex commands
%       for the top, middle, and bottom horizontal lines. Defaults to
%       '\hline' for all three.
%
%       'extra_hlines' - pass either a vector of numeric or logical indices
%       to add extra horizontal lines in the table. They will be added
%       after the rows given by the numeric indices or after rows for which
%       the logical index is true. If given as a logical index, it is
%       checked that it is the same length as the number of rows in the
%       table. lines{2} is used for the horizontal line command.
%
%       'environment' - lets you change the float type. Default is 'table',
%       you could change this to e.g. 'sidewaystable' to put the table
%       sideways on the page (that particular example requires the
%       "rotating" package is loaded in your .tex file).
%
%       'center' - boolean, adds the \centering command to the table
%       environment to ensure the table is centered.

E = JLLErrors;

p = inputParser;
p.addParameter('file', '', @ischar);
p.addParameter('caption', '');
p.addParameter('label','');
p.addParameter('rownames', {});
p.addParameter('multirow', true);
p.addParameter('colnames', {});
p.addParameter('multicol', true);
p.addParameter('m2l', {});
p.addParameter('insert',false);
p.addParameter('marker','');
p.addParameter('overwrite',false);
p.addParameter('extra_hlines',[]);
p.addParameter('lines', {'\hline', '\hline', '\hline'});
p.addParameter('environment', 'table');
p.addParameter('center', false);

p.parse(varargin{:});
pout = p.Results;

file_out = pout.file;
caption = pout.caption;
label = pout.label;
rownames = pout.rownames;
use_multirow = pout.multirow;
colnames = pout.colnames;
use_multicol = pout.multicol;
m2l_opts = pout.m2l;
do_insert = pout.insert;
insert_mark = pout.marker;
overwrite = pout.overwrite;
extra_hlines = pout.extra_hlines;
hlines = pout.lines;
environment = pout.environment;
do_center = pout.center;

if istable(T) && isempty(rownames)
    rownames = T.Properties.RowNames;
end
if istable(T) && isempty(colnames)
    colnames = T.Properties.VariableNames;
end

% Handle any % signs in the row/column names so that they aren't lost in
% the various sprintf calls.
rownames = strrep(rownames, '%', '%%');
colnames = strrep(colnames, '%', '%%');

% Check the input variables are of the right type
if ~ischar(T) && ~isnumeric(T) && ~istable(T)
    E.badinput('T must be a string, array, or table.')
end

if ~ischar(file_out)
    E.badinput('The parameter "file_out" must be a string')
end

if ~ischar(caption)
    E.badinput('The parameter "caption" must be a string')
end

if ~ischar(label)
    E.badinput('The parameter "label" must be a string')
end

if ~iscellstr(rownames)
    E.badinput('The parameter "rownames" must be a cell array of strings')
end

if ~iscellstr(colnames)
    E.badinput('The parameter "colnames" must be a cell array of strings')
end

if ~iscell(m2l_opts)
    E.badinput('The parameter "m2l" must be a cell array')
end

if ~isscalar(do_insert) || (~isnumeric(do_insert) && ~islogical(do_insert))
    E.badinput('The parameter "insert" must be a scalar number or boolean value')
end

if ~ischar(insert_mark)
    E.badinput('The parameter "marker" must be a string')
end

if ~isscalar(overwrite) || (~isnumeric(overwrite) && ~islogical(overwrite))
    E.badinput('The parameter "overwrite" must be a scalar number or boolean value')
end

if ~iscellstr(hlines) || numel(hlines) ~= 3
    E.badinput('The parameter "lines" must be a 3 element cell array of strings')
end

% Finally check the interrelationships
if do_insert && isempty(insert_mark)
    E.badinput('To use ''insert'' == true, you must specify a marker string')
end


%%%%% MAIN FUNCTION %%%%%
is_rownames = ~isempty(rownames) > 0;
if isrow(rownames)
    % Allow the user to give a single column of rownames as either a row or
    % column vector
    rownames = rownames';
end

if ischar(T)
    latex_table_body = T;
elseif isnumeric(T) || istable(T)
    if istable(T)
        T = table2array(T);
    end
    latex_table_body = mat2latex(T, m2l_opts{:});
end

n_table_rows = numel(strfind(latex_table_body, '\\'));
if is_rownames && size(rownames,1) ~= n_table_rows
    % It's possible that a string passed in as T might not have a \\ at the
    % end of the last line. This is a kludgy solution, but I don't use the
    % character input very ofter
    if size(rownames,1) == n_table_rows + 1
        warning('rownames may have one too many rows unless you passed T as a character array and did not include a \\ after the last row')
    else
        E.badinput('If giving row names, the cell array must have the same number of rows as the table (if giving only a single column of row names, it may be a row vector with the same number of columns as the table has rows)')
    end
end

% Check the extra hlines. If it's not logical, convert to a logical array
if isempty(extra_hlines)
    extra_hlines = false(n_table_rows + 1,1); % The extra +1 accounts for if T is passed as a character array and is missing it's final \\. Having extra will not hurt.
elseif ~islogical(extra_hlines)
    tmp = false(n_table_rows + 1, 1);
    tmp(extra_hlines) = true;
    extra_hlines = tmp;
end

if numel(extra_hlines) ~= n_table_rows && numel(extra_hlines) ~= (n_table_rows + 1)
    E.badinput('If giving extra hlines as a logical array, it must have the same number of elements as there are rows in the table')
end

latex_table_body = strsplit(latex_table_body, '\n');
latex_table_body = latex_table_body(~iscellcontents(latex_table_body,'isempty'));

% If no row names, we don't need to reserve a column for them in the
% header, nor print them (obviously)
[rownames, n_columns_rownames] = format_row_names(rownames, use_multirow);

n_columns = numel(strfind(latex_table_body{1},'&')) + 1 + n_columns_rownames;

% Allow the column names to include one for the row names or not.
if ~isempty(colnames)
    if size(colnames,2) == n_columns 
        extra_ands = '';
    elseif size(colnames,2) == n_columns - n_columns_rownames
        extra_ands = repmat('& ', 1, n_columns_rownames);
    else
        E.badinput('Wrong number of column names: if given, the parameter ''colnames'' must include an entry for each data column OR an entry for every data column and row name column')
    end
    
    header = format_header(colnames, use_multicol, extra_ands);
    % Add the middle dividing line before the last newline
    header = regexprep(header, '\\n $', strrep([hlines{2},'\n'],'\','\\'));
    tabular_body = sprintf('%s %s\n ', hlines{1}, header);
else
    tabular_body = [hlines{1}, ' '];
end

for a=1:numel(latex_table_body)
    if a == numel(latex_table_body)
        sep = sprintf('%s \\n ', hlines{3});
    elseif extra_hlines(a)
        sep = sprintf('%s \\n ', hlines{2});
    else
        sep = '\n ';
    end
    
    if is_rownames
        tabular_body = [tabular_body, rownames{a}, ' & ', latex_table_body{a}, sep];
    else
        tabular_body = [tabular_body, latex_table_body{a}, sep];
    end
end
tabular_body = sprintf(tex_in_printf(tabular_body));

% Now figure out how tabular should be formatted
if is_rownames
    fmt_str = sprintf('%s%s',repmat('l',1,n_columns_rownames), repmat('r', 1, n_columns-n_columns_rownames));
else
    fmt_str = sprintf('%s', repmat('r', 1, n_columns));
end

% And whether caption and label should be included
if ~isempty(caption)
    cap_str = sprintf('\\caption{%s}\n ', caption);
else
    cap_str = '';
end

if ~isempty(label)
    label_str = sprintf('\\label{%s}\n ', label);
else
    label_str = '';
end

% Put it all together
if do_center
    center_line = ' \centering\n';
else
    center_line = ' ';
end
table_str = [' \begin{', environment,'}\n',...
             center_line,...
             ' \begin{tabular}{%1$s}\n'...
             ' %2$s\n'...
             ' \end{tabular}\n'...
             ' %3$s%4$s'...
             '\end{', environment,'}\n'];
table_str = tex_in_printf(table_str);
    
if isempty(file_out)
    varargout{1} = sprintf(table_str, fmt_str, tabular_body, cap_str, label_str);
elseif ~do_insert    
    if exist(file_out, 'file') && ~overwrite
        user_ans = ask_yn(sprintf('File %s exists. Overwrite?', file_out));
        if ~user_ans
            return
        end
    end
    fid = fopen(file_out,'w');
    if fid < 0
        E.callError('io_error', 'Could not open %s for writing', file_out);
    end
    fprintf(fid, table_str, fmt_str, tabular_body, cap_str, label_str);
    fclose(fid);
else
    insert_into_file(sprintf(table_str, fmt_str, tabular_body, cap_str, label_str), file_out, insert_mark);
end

end

function insert_into_file(table_str, file_name, mark)
E = JLLErrors;
table_str = strrep(table_str, '%', '%%');

bckp_name = sprintf('%s.bckp', file_name);
out_name = sprintf('%s.matlab', file_name);

if ~exist(file_name, 'file')
    E.filenotfound(file_name);
elseif exist(file_name, 'dir')
    E.callError('file_is_dir', 'Given file (%s) is actually a directory', msg);
end

[stat, msg] = copyfile(file_name, bckp_name);
if ~stat
    E.callError('could_not_copy','Could not make backup copy: %s', msg);
end

[curr_fid, msg] = fopen(file_name, 'r');
if curr_fid < 0
    E.callError('could_not_open','Could not open %s for reading: %s', file_name, msg);
end

[new_fid, msg] = fopen(out_name, 'w');
if new_fid < 0
    E.callError('could_not_open','Could not open %s for writing: %s', file_out, msg);
end

found_mark = false;
in_mark = false;
mark_regex_start = sprintf('^%%BEGIN\\s+%s',mark);
mark_regex_end = sprintf('^%%END\\s+%s',mark);
mark_regex = mark_regex_start;
tline = fgets(curr_fid);
while ischar(tline)
    nline = strrep(tline, '%', '%%');
    if ~in_mark
        fprintf(new_fid, tex_in_printf(nline));
    end
    if ~isempty(regexp(tline, mark_regex, 'once'))
        if ~in_mark
            if found_mark
                warning('Mark %s multiply defined', mark);
            end
            found_mark = true;
            in_mark = true;
            mark_regex = mark_regex_end;
        elseif in_mark
            fprintf(new_fid, tex_in_printf(table_str));
            fprintf(new_fid, tex_in_printf(nline));
            in_mark = false;
            mark_regex = mark_regex_start;
        end 
    end
    
    tline = fgets(curr_fid);
end
fclose(curr_fid);
fclose(new_fid);

if in_mark
    E.callError('bad_mark', 'Could not find ending mark line for %s', mark);
elseif ~found_mark
    E.callError('no_mark', 'Could not find mark %s', mark);
end

movefile(out_name, file_name);

end

function header = format_header(colnames_in, use_multicol, extra_ands)
% Take the column names cell array and format it for inclusion in the latex
% table. This means joining columns with the \& separator and lines with
% the \\ separator. If USE_MULTICOL is true, then adjacent cells containing
% the same string will be merged into a single multicol function. Currently
% it will always center the header within each column and does not include
% lines to either side.
header = '';
if use_multicol
    multicol = 'col';
else
    multicol = 'no';
end

for i_row = 1:size(colnames_in,1)
    header = [header, extra_ands, format_header_row(colnames_in(i_row,:), multicol, true), ' \\ \n '];
end
end

function [row_names, n_columns] = format_row_names(row_names_in, use_multirow)
% If row_names_in is empty, there's nothing to do (and no columns to
% reserve)
if isempty(row_names_in)
    row_names = row_names_in;
    n_columns = 0;
    return
end

if use_multirow
    multirow = 'row';
else
    multirow = 'no';
end

% Transpose the row names so that we can reuse the code in
% format_header_row, which is written for the header and combines identical
% entries adjacent along rows.
row_names_in = row_names_in';
for i_row = 1:size(row_names_in,1)
    row_names_in(i_row,:) = format_header_row(row_names_in(i_row,:), multirow, false);
end

% Now we transpose back and join each true row of the row names into a
% single string connected by &'s
row_names_in = row_names_in';
row_names = cell(size(row_names_in,1),1);
for i_row = 1:numel(row_names)
    row_names{i_row} = strjoin(row_names_in(i_row,:), ' & ');
end

n_columns = size(row_names_in,2);

end

function row = format_header_row(row, use_multi, do_join)
E = JLLErrors;

if ~iscellstr(row)
    E.badinput('ROW must be a cell array of char arrays');
end

allowed_use_multi = {'no','col','row'};
if ~ismember(use_multi, allowed_use_multi)
    E.badinput('USE_MULTI must be one of: %s', strjoin(allowed_use_multi, ', '));
end

i_col = 1;
keep_columns = true(size(row));
while ~strcmpi(use_multi,'no') && i_col <= numel(row)
    i_col_same = 1;
    
    while i_col + i_col_same <= numel(row) && strcmp(row{i_col}, row{i_col + i_col_same})
        % Mark the same columns for removal in multicolumn mode. For
        % multirow, each row needs the same number of &'s, so don't remove
        % them, but do replace them repeated entries with empty strings.
        if strcmpi(use_multi, 'col')
            keep_columns(i_col + i_col_same) = false;
        elseif strcmpi(use_multi, 'row')
            row{i_col + i_col_same} = '';
        end
        i_col_same = i_col_same + 1;
    end
    
    % Only replace the current string with a multicolumn if there are
    % actually adjacent identical headers
    if i_col_same > 1
        if strcmpi(use_multi, 'col')
            row{i_col} = sprintf('\\multicolumn{%d}{c}{%s}', i_col_same, row{i_col});
        elseif strcmpi(use_multi, 'row')
            row{i_col} = sprintf('\\multirow{%d}{*}{%s}', i_col_same, row{i_col});
        end
    end
    
    i_col = i_col + i_col_same;
end

if do_join
    row = strjoin(row(keep_columns), ' & ');
end
end
