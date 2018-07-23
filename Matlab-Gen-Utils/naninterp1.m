function [ q ] = naninterp1( xv, v, xq )
%NANINTERP1 Interpolate in 1D ignoring NaNs in XV
%   Q = NANINTERP1( XV, V, XQ ) Interpolates values V defined at
%   coordinates XV to coordinates XQ, first removing NaNs in XV and their
%   corresponding values in V. If there are not at least 2 non-NaN values
%   of XV, then Q will be returned as all NaNs, since interpolation
%   requires at least two points in a given dimension.

E = JLLErrors;

if ~isvector(xv) || ~isvector(v) || ~isvector(xq)
    E.badinput('All inputs must be vectors')
end

notnans = ~isnan(xv);
if sum(notnans) < 2
    q = nan(size(xq));
    return
end

q = interp1(xv(notnans), v(notnans), xq);

end

