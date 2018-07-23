function [ m, sm ] = lsqfitnmorg( x, y )
%LSQFITNMORG Do a least-squares fit of data, passing through the origin
%   [ M, SM ] = LSQFITNMORG( X, Y ) Fits a line to data defined by X, Y,
%   minimizing the orthogonal distance from the line to each point (x_i,
%   y_i) using a SVD decomposition, and forcing the line to pass through
%   the origin. This is best when X and Y have comparable uncertainties and
%   there is good reason to assume no offset, i.e. Y = 0 when X = 0.
%   Returns the slope, M, and the standard deviation of that slope, SM.


% https://www.mathworks.com/matlabcentral/newsreader/view_thread/172430
[~,~,V] = svd([x(:),y(:)],0);
m = -V(1,2)/V(2,2);

% For R2: https://online.stat.psu.edu/~ajw13/stat501/SpecialTopics/Reg_thru_origin.pdf
% except I don't think that it's valid for this sort of regression, that
% seems to be for y-residual regression

% http://www.statisticshowto.com/find-standard-error-regression-slope/ and
% http://stattrek.com/regression/slope-confidence-interval.aspx?Tutorial=AP
yhat = m .* x(:);
sm = sqrt( sum( (y(:) - yhat).^2 ) / (numel(x)-2) ) ./ sqrt( sum(x(:) - mean(x)).^2 );
end

