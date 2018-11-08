function dvec = make_datevec(start_dates, end_dates)
%MAKE_DATEVEC Construct a vector of non-contiguous dates
%
%   DVEC = MAKE_DATEVEC( START_DATES, END_DATES ) will construct
%   a vector of date numbers for the date range or ranges specified
%   by START_DATES and END_DATES. These may be date numbers, vectors
%   of date numbers, date strings, or cell arrays of date strings.
%
%   Examples:
%
%       MAKE_DATEVEC( '2012-01-01', '2012-01-03' ) returns a vector
%       of date numbers 734869:734871 i.e. 1 Jan 2012 to 3 Jan 2012.
%
%       MAKE_DATEVEC( 734869, 734871 ) yields the same as before.
%
%       MAKE_DATEVEC( {'2012-01-01','2012-06-01','2012-09-01'}, 
%       {'2012-01-03','2012-06-03','2012-09-03'} ) returns a vector
%       containing 1-3 Jan 2012, 1-3 Jun 2012, and 1-3 Sept 2012.
%       In general, to create a vector spanning N non-contiguous
%       time periods, START_DATES and END_DATES must be N long cell
%       array vectors where the start and end dates for the i-th period
%       are START_DATES{i} and END_DATES{i}.
%
%       Alternately, the periods may be specified through date number 
%       vectors, i.e. the previous example is equivalent to:
%       MAKE_DATEVEC( [734869,735021,735113],[734871,735023,735115] ).

start_dnums = validate_date(start_dates);
end_dnums = validate_date(end_dates);

dvec = [];
for i = 1:numel(start_dnums)
    dvec = veccat(dvec, start_dnums(i):end_dnums(i));
end

end

