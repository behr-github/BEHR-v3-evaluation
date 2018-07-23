function [ x_out, y_out ] = average_same_x( x,y )
%average_same_x(x,y) Find replicate x values and averages their y's
%   This function finds all instances of replicate x values and averages
%   the corresponding y values.  It will return a sorted version of x and y
%   with x monotonically increasing.  This is primarily intended to prepare
%   data for interpolation. x and y must both be vectors

if isrow(x); x = x'; end
if isrow(y); y = y'; end

S = [x y];

S = sortrows(S);
u = unique(S(:,1));

x_out = zeros(1,numel(u));
y_out = zeros(1,numel(u));

for a=1:numel(u)
    xx = S(:,1)==u(a);
    x_out(a) = u(a);
    y_out(a) = nanmean(S(xx,2));
end

end

