function xx = quality_flag_filter(bad_flags, qflag_in, all_or_any )

% quality_flag_filter - filters an arbitrary number of data sets based on
% their bitarray quality flags.
%
%   This function takes as arguments 1) a number representing the bits that
%   should NOT be set, 2) a list of the quality flags that correspond the
%   the data in, and 3+) each of the data sets that you want filtered.  
%
%   This function works by comparing all quality flags associated with a
%   particular measurement against the bad_flags input with bitand.  If any
%   of the flags match, then that measurement is rejected.  As an example,
%   consider data with a 4-bit flag.  If bad_flags is set to 2 (i.e. 0010),
%   then any measurement with a quality flag where the second bit (**1*) is
%   1 will be removed.
%
%   This returns a logical matrix that is 1 for cells without the flags
%   set, and 0 otherwise.

E = JLLErrors;

if nargin < 3;
    all_or_any = 'all';
elseif ~any(strcmpi({'all','any'},all_or_any));
    E.badinput('all_or_any must be the string ''all'' or ''any''');
else
    % Makes the switch-case statement effectively case insensitive later
    all_or_any = lower(all_or_any);
end

% Check the qflags passed in; if not in a cell array, make them so
if ~iscell(qflag_in);
    qflags = mat2cell(qflag_in);
else
    qflags = qflag_in;
end

% Prep a logical index variable that will keep track of which measurements
% to remove
xx = true(size(qflags));

% Create a bit mask that will be used to check each flag. It should be the
% same class as the flags (e.g. uint16, single, double, int8).

bitmask = eval(sprintf('%s(0)',class(qflag_in{1})));
for b=1:numel(bad_flags)
    bitmask = bitset(bitmask,bad_flags(b));
end

for a=1:numel(qflags)
    if any(bitand(bitmask,qflags{a}))
        xx(a) = false;
    end
end

end
