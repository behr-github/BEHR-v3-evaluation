function [ avg_slope, avg_int, avg_R2 ] = bootstrapping_analysis( x, y, number_of_trials, varargin )
%bootstrapping_analysis Randomly selects a portion of the data and does a
%fit, repeated n times.
%
%   This requires 3 inputs: x and y values for the data, and the number of
%   times the sampling should be carried out (see below for explanation).
%   10000 is a good first guess value.
%
%   Performs a bootstrapping analysis on the linear regression of data
%   (x,y).  By default, this functions by resampling the data (with
%   replacement) at random until a resample of the same size as the initial
%   data is obtained; a regression fit is then performed.  This process is
%   repeated number_of_trials times, and the resulting slope, intercept,
%   and R^2 values averaged.  
%
%   Additional parameters are:
%       fraction = What fraction of the original sample size the resample
%           size should be.  Defaults to 1; i.e. a sample as big as the
%           original x,y vectors will be created.  If replacement is off,
%           then this should be set to < 1, with 0.5 as a good first guess.
%
%       replacement = Whether or not to allow the algorithm to pick the
%           same point more than once.  Defaults to 1; set to 0 to cause
%           random sampling without replacement.
%
%       fit_type = 'RMA, 'MA' or 'y-resid'; RMA uses a major axis fit that
%           allows for variation in x- and y-values, and allows the scale
%           of the x- and y- axes to be different.  MA uses a simpler
%           regression that also allows both x and y to vary, but requires
%           that their scale be the same (see info from
%           http://www.mbari.org/staff/etp3/regress/index.htm on lsqfitma
%           and lsqfitgm.  'y-resid' uses a common linear regression that
%           only minimized difference in the y values from the fit.  
%
%   This requires the functions lsqfitgm, lsqfitx, and lsqfity from
%   http://www.mbari.org/staff/etp3/regress/index.htm for RMA regression,
%   lsqfitrm from the same site for RM regression, and polyfit_R2 for
%   y-residual fits.  polyfit_R2 is an extension of the built-in polyfit
%   function that also calcuates R^2 values.

p = inputParser;
p.addRequired('x',@isvector);
p.addRequired('y',@isvector);
p.addRequired('number_of_trials',@isscalar);
p.addParamValue('fraction',1,@isscalar);
p.addParamValue('replacement',1,@isscalar);
p.addParamValue('fit_type','RMA',@isstr);

p.parse(x,y,number_of_trials,varargin{:});
pout = p.Results;

x = pout.x;
y = pout.y;
number_of_trials = pout.number_of_trials;
percent_of_data = pout.fraction;
replacement = pout.replacement;
fit_type = pout.fit_type;

if isrow(x); x = x'; end
if isrow(y); y = y'; end
if numel(x) ~= numel(y); error('bootstrap:xy_unequal','x and y must have equal lengths'); end
if ~any(strcmpi(fit_type,{'RMA','MA','y-resid'}));
    error('bootstrap:fit_type','Fit type must be ''RMA'' (reduced major axis), ''MA'' (major axis) or ''y-resid'' (typical linear regression, y-residuals only \n See help for details on selection.');
end

nans = isnan(x) | isnan(y);
x(nans) = []; y(nans) = [];
M = [x,y];
count = ceil(percent_of_data*numel(x));

slopes = zeros(number_of_trials,1);
ints = zeros(number_of_trials,1);
Rsquareds = zeros(number_of_trials,1);

for a=1:number_of_trials
    x_set = zeros(count,1); y_set = zeros(count,1);
    M_set = M;
    for b=1:count
        ri = randi([1, size(M_set,1)],1);
        x_set(b) = M_set(ri,1); y_set(b) = M_set(ri,2);
        if replacement < 1; M_set(ri,:) = []; end
    end
    
    switch fit_type
        case 'RMA'
            [m,b,r] = lsqfitgm(x_set,y_set);
            r=r^2;
        case 'MA'
            [m,b,r] = lsqfitma(x_set,y_set);
            r = r^2;
        case 'y_resid'
            [P,r] = polyfit_R2(x_set,y_set,1);
            m = P(1);
            b = P(2);
    end
    slopes(a) = m;
    ints(a) = b;
    Rsquareds(a) = r;
end

avg_slope = mean(slopes);
avg_int = mean(ints);
avg_R2 = mean(Rsquareds);
end

