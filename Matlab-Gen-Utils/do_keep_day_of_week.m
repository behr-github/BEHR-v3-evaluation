function [ keep_bool ] = do_keep_day_of_week(curr_date, days_of_week)
%DO_KEEP_DAY_OF_WEEK Should the current date be kept based on its day-of-week?
%   KEEP_BOOL = DO_KEEP_DAY_OF_WEEK( CURR_DATE, DAYS_OF_WEEK ) Returns true
%   if the weekday of CURR_DATE is one of those specified in DAYS_OF_WEEK,
%   false otherwise. CURR_DATE may be any date number or date string
%   understood by Matlab. DAYS_OF_WEEK must be a character array with
%   letters representing days of the week to keep:
%       U = Sunday
%       M = Monday
%       T = Tuesday
%       W = Wednesday
%       R = Thursday
%       F = Friday
%       S = Saturday

E = JLLErrors;

validate_date(curr_date);
single_letter_abbrevs = {'U','M','T','W','R','F','S'};
if ~check_days_of_week(days_of_week, single_letter_abbrevs)
    E.badinput('DAYS_OF_WEEK must be a character array consisting of only the letters U, M, T, W, R, F, or S in some combination');
end

numeric_days_of_week = map_days_to_number(days_of_week, single_letter_abbrevs);
keep_bool = any(weekday(curr_date) == numeric_days_of_week);

end

function days_as_numbers = map_days_to_number(days_of_week_letters, single_letter_abbrevs)
% This maps single letter abbreviations of days of the week to the numbers
% returned by weekday(), which gives Sunday = 1 and Saturday = 7.
% days_as_numbers will be the numeric equivalent of days_of_week_letters.

tmp_bool = false(size(single_letter_abbrevs));

for i_day = 1:numel(single_letter_abbrevs)
    tmp_bool(i_day) = ismember(single_letter_abbrevs{i_day}, days_of_week_letters);
end

days_as_numbers = find(tmp_bool);

end

function chk = check_days_of_week(days_of_week, allowed_abbrevs)
chk = false;

if ~ischar(days_of_week)
    return
end

for i_day = 1:numel(days_of_week)
    if ~ismember(days_of_week(i_day), allowed_abbrevs)
        return
    end
end

chk = true;
end