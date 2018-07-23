function varargout = square_subplot_dims(nplots, varargin)
%SQUARE_SUBPLOT_DIMS Calculate best subplot dimensions
%   SP_DIMS = SQUARE_SUBPLOT_DIMS( NPLOTS ) Calculates the best subplot
%   dimensions for the number of plots, NPLOTS, trying to keep the layout
%   as square as possible. Prefers wide over tall (more plots horizontally
%   than vertically) if an even distribution cannot be obtained.
%
%   SP_DIMS = SQUARE_SUBPLOT_DIMS( ___, 'tall' ) will prefer tall over wide
%   instead.
%
%   SP_DIMS = SQUARE_SUBPLOT_DIMS( ___, 'exact' ) will try to get the exact
%   right number of subplots, even if it makes the figure more rectangular.
%   However, it will only go so far, and may still give dimensions that do
%   not exactly equal the number of plots if to do so would require a very
%   rectangular shape.
%
%   [ H, W ] = SQUARE_SUBPLOT_DIMS( ___ ) returns the subplot dimensions
%   separately instead of as a single vector. Returns height, then width.

p = advInputParser;
p.addFlag('tall');
p.addFlag('exact');
p.parse(varargin{:});
pout = p.AdvResults;

prefer_tall = pout.tall;
exact = pout.exact;

sz(1) = ceil(sqrt(nplots));
sz(2) = ceil(nplots / sz(1));

% If we don't have exactly the right number of plots, try another option.
% What if we had used floor() instead?
if prod(sz) ~= nplots && exact
    alt = ceil(nplots / (sz(1) - 1));
    % If this doesn't get us exactly the right number of plots, go back to
    % the first way.
    if (sz(1) - 1)*alt == nplots
        sz = [sz(1) - 1, alt];
    end
end

if prefer_tall
    sp_dims = [max(sz), min(sz)];
else
    sp_dims = [min(sz), max(sz)];
end

if nargout < 2
    varargout = {sp_dims};
else
    varargout = num2cell(sp_dims);

end

