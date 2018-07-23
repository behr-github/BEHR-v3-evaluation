function dvec = make_datevec(start_dates, end_dates)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

start_dnums = validate_date(start_dates, 'cell');
end_dnums = validate_date(end_dates, 'cell');

dvec = [];
for i = 1:numel(start_dnums)
    dvec = veccat(dvec, start_dnums(i):end_dnums(i));
end

end

