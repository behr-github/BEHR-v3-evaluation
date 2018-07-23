function [ r ] = randdist( pdffxn, xrange, sz )
%R = RANDDIST( PDFFXN, XRANGE ) Return a random number, R, from the defined PDF
%   This function will return a random number, R, drawn from the
%   distribution defined by the function handle PDFFXN.  This should be a
%   function that given a single vector input x returns a vector of the
%   same size as x. You must also specify the range of valid values using
%   the two element vector XRANGE, as this function does not currently look
%   for where the PDF drops to approximately 0. You can also use this to
%   impose a soft limit on the values the random number can take on.
%
%   RANDDIST( PDFFXN, XRANGE, SZ ) will return R as an array with with size
%   defined by the vector SZ.
%
%   Examples:
%   1) Draw a number 1000 times from a normal distribution centered on
%   0.5:
%       
%       t = nan([1 1000]);
%       f = @(x) exp(-(x-0.5).^2 ./ 0.05);
%       for i=1:1000
%           t(i) = randdist(f, [0 1]);
%       end
%       figure; 
%       histogram(t,20,'normalization','pdf');
%       x=0:0.01:1;
%       line(x, f(x), 'color', 'r')
%
%   2) Draw a 10,000 element vector from a flat-topped Gaussian
%   distribution centered on 10:
%
%       g = @(x) exp( -(x-10).^12 );
%       t = randdist(g, [8 12], [1 10000]);
%       figure;
%       histogram(t,50,'normalization','pdf');
%       x = 8:0.01:12;
%       line(x, g(x), 'color', 'r');

E = JLLErrors;

if ~exist('sz','var')
    sz = 1;
end

x = linspace(min(xrange), max(xrange), 500);
cdfvec = cumsum(pdffxn(x));

if ~isequal(size(cdfvec), size(x))
    E.badinput('The PDF function must be such that given a 1-by-n vector, it also returns a 1-by-n vector')
end

r = nan(sz);
for i1=1:numel(r)
    r0 = rand;
    cdffxn = cdfvec ./ max(cdfvec) - r0;
    
    % Find between which two points the CDF - rand crosses zero.
    for a=2:numel(cdffxn)
        if sign(cdffxn(a)) ~= sign(cdffxn(a-1))
            k=a;
            % Now compute the exact point at which it crosses
            m = (cdffxn(k) - cdffxn(k-1))/(x(k) - x(k-1));
            b = cdffxn(k) - m*x(k);
            r(i1) = -b/m;
            break
        elseif a==numel(cdffxn)
%            E.callError('no_zero','Found no zero crossing point. CDF not bounded on 0 to 1?');
            % Assuming that the zero crossing is somewhere between 0 and
            % the first value of CDFFXN. Further assuming that the spacing
            % in x is relatively constant
            m = cdffxn(1) / (x(2)-x(1));
            b = cdffxn(1) - m*x(1);
            r(i1) = -b/m;
        end
    end
    
    
end
end

