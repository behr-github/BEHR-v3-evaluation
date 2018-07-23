function [pandora_dnum] = pandora_gmt_to_datenum(pandora_gmt)
%PANDORA_GMT_TO_DATENUM Convert Pandora time values to date numbers
%   PANDORA_DNUM = PANDORA_GMT_TO_DATENUM( PANDORA_GMT ) Converts the
%   Pandora GMT timestamp PANDORA_GMT to a date number, PANDORA_DNUM.
%   PANDORA_GMT may be a numeric array, cell array of character arrays, or
%   (if only a single date) a character array. In any case, it is expected
%   to be formatted as yyyymmddHHMMSS.

E = JLLErrors;

if isnumeric(pandora_gmt)
    tmp = cell(size(pandora_gmt));
    for i_date = 1:numel(pandora_gmt)
        tmp{i_date} = num2str(pandora_gmt(i_date));
    end
    pandora_gmt = tmp;
elseif ischar(pandora_gmt)
    pandora_gmt = {pandora_gmt};
elseif ~iscellstr(pandora_gmt)
    E.badinput('PANDORA_GMT must be a numeric array or cell array of char arrays');
end

pandora_dnum = nan(size(pandora_gmt));
for i_date = 1:numel(pandora_gmt)
    pandora_dnum(i_date) = datenum(pandora_gmt{i_date}, 'yyyymmddHHMMSS');
end

end

