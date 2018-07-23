function cbdatetick(cb, varargin)
%CBDATETICK Label a colorbar with dates
%   CBDATETICK( CB ) Changes the tick labels on the colorbar referenced by
%   handle CB to dates in YYYY-MM-DD format, assuming that the color values
%   are date numbers.
%
%   CBDATETICK( CB, DATEFMT ) Uses DATEFMT as the second argument to
%   DATESTR() to convert the date number represented by each tick mark to a
%   human-readable date string.

p = advInputParser;
p.addOptional('datefmt', 'yyyy-mm-dd', @(x) ischar(x) || isstring(x))

p.parse(varargin{:});
pout = p.Results;
date_format = pout.datefmt;

tick_values = cb.Ticks;
tick_labels = cell(size(cb.Ticks));

for i_tick = 1:numel(tick_values)
    tick_labels{i_tick} = datestr(tick_values(i_tick), date_format);
end

cb.TickLabels = tick_labels;

end

