function [ dbout, varargout ] = match_arrays2db( db, varargin )
%match_arrays2db Match values from arrays and db structure
%   My methods that run both boundary layer and spiral validation output
%   results in two forms: 1) cell arrays that contain the main lat, lon,
%   and column data, 2) a structure called "db" that contains a bunch of
%   addition and debugging information.  Each cell in the arrays matches to
%   one top level index of db, but matching individual measurements is not
%   easily accomplished with a few built in Matlab functions.  This
%   function will do that for you.
%
%   As input this takes the db structure plus whatever of the output cell
%   arrays you want matched. It will return variables in the same order.
%   The cell arrays will be converted to column vectors that contain all
%   the values in the cell arrays extracted into the vectors.  The db
%   structure will be returned as a structure, but with only one top-level
%   index.  Each field will be in one of two formats:
%       1) A cell array where cell{i} contains all the values for that
%       field that correspond to entry(i) in the vectors returned.
%       2) A column vector where entry(i) is the mean or median of the
%       values that correspond to entry(i) in the main output vectors.
%   By default, option 1 is used.  To use option 2 instead, pass 'mean' or
%   'median' as a string as the second or later argument.  An exception is
%   a lat/lon corner field (i.e. any field name containing the substring
%   'corn' and 'lat' or 'lon').  These will not be averaged and will always
%   be returned in a n x 4 matrix.

argin = varargin;

% Default to outputing db fields as cell arrays, but if the user passes a
% string, use that method instead.
method = 'cell';
for a=1:numel(argin)
    if ischar(argin{a});
        method = argin{a};
        argin(a) = [];
    end
end

% To initialize the variables, assume that on average there are two entries
% per cell
n = 2*numel(db);
fields = fieldnames(db);
dbout = struct;
argout = repmat({zeros(n,1)},1,numel(argin));
for a=1:numel(argin)
    if iscell(argin{a}{1})
        argout{a} = cell(n,1);
    end
end

if strcmpi(method,'cell'); dummy = {cell(n,1)};
elseif any(strcmpi(method,{'mean','median'})); dummy = zeros(n,1);
else error('match_arrays2db:method','Method string must be ''cell'', ''median'', or ''mean''');
end
for a=1:numel(fields)
    if cornerfield(fields{a});
        dbout.(fields{a}) = {cell(4,n)};
    else
        dbout.(fields{a}) = dummy;
    end
end

E=0;
% Loop through each top level entry in db and each cell in the arrays
for a=1:numel(db)
    % Loop through each entry within those variables
    for b=1:numel(db(a).(fields{1}))
        E=E+1;
        % Loop through each field and array.  If it is a lat/lon corner
        % field, never average it 
        for c=1:numel(fields)
            if cornerfield(fields{c})
                dbout.(fields{c}){E} = db(a).(fields{c}){b};
            else
                switch method
                    case 'cell'
                        dbout.(fields{c}){E} = db(a).(fields{c}){b};
                    case 'mean'
                        dbout.(fields{c})(E) = nanmean(db(a).(fields{c}){b});
                    case 'median'
                        dbout.(fields{c})(E) = nanmedian(db(a).(fields{c}){b});
                end
            end
        end
        for c=1:numel(argin);
            argout{c}(E) = argin{c}{a}(b);
        end
    end
end

% Clean up the output, removing the extra values.
if E < n;
    for a=1:numel(fields);
        dbout.(fields{a}) = dbout.(fields{a})(1:E);
    end
    for a=1:numel(argout)
        argout{a} = argout{a}(1:E);
    end
end

varargout = argout;


function test = cornerfield(field)
t1 = ~isempty(regexpi(field,'lat'));
t2 = ~isempty(regexpi(field,'lon'));
t3 = ~isempty(regexpi(field,'corn'));

% Return true if the field name contains either 'lat' or 'lon' as well as
% 'corn' (test for pixel corner field names)
test = (t1 || t2) && t3;

