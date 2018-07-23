function [ p_val ] = p_val_slope( x, y, regression_type, one_tailed_bool )
%p_val_slope Calculate the p-value for a linear regression
%   The p-value represents the probability that the slope is actually equal
%   to 0. To calculate this, I will follow the method in:
%
%       http://www.real-statistics.com/regression/hypothesis-testing-significance-regression-line-slope/
%   
%   which calculates a test t statistic and then compares that to a
%   t-distribution.  This assumes that your data can be adequately
%   represented by a t-distribution (i.e. is normal and unskewed). 
%
%   This function requires 3 inputs, and can function in one of two ways:
%       I) Given two vectors (x & y) and a regression type as a string
%       (y-resid, x-resid, majoraxis, or RMA), this will calculate the
%       p-value based on a fit for that data set.
%
%       II) Given a slope, standard deviation, and a number of points, this
%       will use those directly to calculate the p-value. 
%
%   The final argument is optional and should be 1 or 0. If 1, this will
%   treat the distribution as one-tailed, i.e. the data is only on one side
%   of 0. This is actually a less strict criterion than the two-tailed test
%   because it will allocate all of the probability that the slope is
%   different from 0 to one side of the distribution, instead of both
%   sides. See:
%
%       http://blog.minitab.com/blog/adventures-in-statistics/understanding-hypothesis-tests%3A-significance-levels-alpha-and-p-values-in-statistics
%
%   This defaults to 0, i.e. this will report a two-tailed p-value. You
%   should only use a one-tailed test if you are confident that the slope
%   should be positive or negative (and not it could be either). An example
%   would be comparing satellite and aircraft VCDs - we would expect there
%   to be a positive correlation.  If in doubt, use the two-tailed test.
%
%   Other resources:
%       http://blog.minitab.com/blog/adventures-in-statistics/how-to-correctly-interpret-p-values
%       http://stattrek.com/regression/slope-test.aspx?Tutorial=AP
%       http://www.ats.ucla.edu/stat/mult_pkg/faq/general/tail_tests.htm
%
%   Note that this function is only
%
%   Josh Laughner <joshlaugh5@gmail.com> 23 Apr 2015

E = JLLErrors;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT VALIDATION %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    one_tailed_bool = 0;
end

% First, check for the input mode where x & y are the slope and std. error
% rather than vectors of data
if isscalar(x) && isscalar(y)
    if ~isscalar(regression_type)
        E.badinput('The third input must be a scalar NUMBER OF POINTS (not degrees of freedom) if inputting slope/std. error');
    end
    
    % We'll rename the variables to be clearer
    slope = x;
    stderr = y;
    numpts = regression_type;
    
    % Use this to determine which way to calculate the p-value
    slope_in_bool = true;
elseif isvector(x) && isvector(y)
    if ~all(size(x) == size(y))
        E.badinput('x & y must be the same size, if they are vectors');
    elseif ~ischar(regression_type) || ~any(strcmpi(regression_type,{'y-resid','x-resid','majoraxis','rma'}))
        E.badinput('regression_type must be a string: y=resid, x-resid-, majoraxis, or rma');
    end
    
    % Use this to determine which way to calculate the p-value
    slope_in_bool = false;
        
else
    E.badinput('The first two inputs must either both be vectors of data or a scalar value of slope and std. error of slope');
end

if ~isscalar(one_tailed_bool)
    E.badinput('one_tailed_bool must be a scalar (0 or 1)');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% CALCULATE P-VALUE %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If passed data, use the calc_fit_line function to get the slope and
% standard error
if ~slope_in_bool
    [~,~,~,S] = calc_fit_line(x,y,'regression',regression_type);
    slope = S.P(1);
    stderr = S.StdDevM;
    numpts = sum(~isnan(x) & ~isnan(y));
end

t = slope / stderr;
if ~isscalar(t)
    E.badvartype(t,'scalar');
end

% If the t-value is negative (slope is negative) then we should take the
% p-value as the cumulative probability up to the value of t. If t > 0,
% then we want the p value as the cumulative probability not yet accounted
% for.

p_val = tcdf(t,numpts);

if t > 0
    p_val = 1 - p_val;
end

% If we're doing a two-tailed test, then the actual area under the curve we
% want needs to be multiplied by 2 to account for the other side.
if ~one_tailed_bool
    p_val = 2 * p_val;
end

end

