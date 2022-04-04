% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 4/2/22
%
% This script creates the Vandermonde matrix to interpolate cos(x) on the
% interval [0, pi/2] by a fourth order polynomial. We will then evaluate
% this polynomial at random points in this interval to see if the error
% obeys a derived upper bound

p0 = @(x) 1;
p1 = @(x) x;
p2 = @(x) x^2;
p3 = @(x) x^3;
p4 = @(x) x^4;

% points at which to interpolate
x = [0, pi/8, pi/4, 3*pi/8, pi/2];
y = cos(x)';

% vandermonde matrix
A = [arrayfun(p0, x); arrayfun(p1, x); arrayfun(p2, x); arrayfun(p3, x); arrayfun(p4, x)]';

% solve for the coefficients
c = A\y;

% interpolating polynomial
p = @(x) c(1)*p0(x)+c(2)*p1(x)+c(3)*p2(x)+c(4)*p3(x)+c(5)*p4(x);

