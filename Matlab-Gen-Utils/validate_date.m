function [ varargout ] = validate_date( date_in )
%VALIDATE_DATE Check if a date can be understood
%
%   VALIDATE_DATE( DATE_IN ) Checks in DATE_IN can be understood by Matlab.
%   If it is numeric, it checks that it can be converted to a date string;
%   if a string, checks that it can automatically be converted to a date
%   number. If a cell array, checks that each cell's contents follow the
%   above rules.
%
%   DATE_OUT = VALIDATE_DATE( DATE_IN ) Returns the validated date(s) as a
%   date number or, if DATE_IN is a cell array, an array of date numbers
%   the same shape as DATE_IN.

if iscell(date_in)
    date_out = nan(size(date_in));
    for a=1:numel(date_in)
        date_out(a) = validate_single_date(date_in{a});
    end
else
    date_out = validate_single_date(date_in);
end

if nargout > 0
    varargout{1} = date_out;
end

end

function date_out = validate_single_date(date_in)
E = JLLErrors;
if ischar(date_in)
    try 
        date_out = datenum(date_in);
    catch err
        if strcmp(err.identifier, 'MATLAB:datenum:ConvertDateString')
            E.baddate(date_in);
        else
            rethrow(err)
        end
    end
elseif isnumeric(date_in)
    try
        datestr(date_in);
    catch err
        if strcmp(err.identifier, 'MATLAB:datestr:ConvertDateNumber')
            E.baddate('Numerical date (%g) cannot be understood by datestr - value is invalid', date_in)
        else
            rethrow(err)
        end
    end
    date_out = date_in;
else
    E.baddate(date_in);
end
end