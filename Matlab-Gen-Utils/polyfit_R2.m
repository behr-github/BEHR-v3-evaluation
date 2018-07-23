function [ P, R2 ] = polyfit_R2( x, y, n )
%polyfit_R2(x,y,n): Fits x & y to a polynomial of degree n and returns an R-squared value.

P = polyfit(x,y,n);

yfit = polyval(P,x);
yresid = y - yfit;
SSresid = sum(yresid.^2);
SStotal = (numel(y)-1) * var(y);
R2 = 1 - SSresid/SStotal;
end

