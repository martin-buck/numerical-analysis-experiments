% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 4/2/22
%
% This script constructs the Newton and Lagrange basis polynomials to
% estimate the drag coefficient d(v) which is a function of velocity for a
% given set of input-output pairs, and uses the approximating polynomial to
% estimate the drag at an interpolating point

% The data to be interpolated
x = [0, 50, 75, 100, 125];
y = [.5, .5, .4, .28, .23];

% The Lagrange polynomials
l1 = @(t) ((t-x(2))*(t-x(3))*(t-x(4))*(t-x(5)))/((x(1)-x(2))*(x(1)-x(3))*(x(1)-x(4))*(x(1)-x(5)));
l2 = @(t) ((t-x(1))*(t-x(3))*(t-x(4))*(t-x(5)))/((x(2)-x(1))*(x(2)-x(3))*(x(2)-x(4))*(x(2)-x(5)));
l3 = @(t) ((t-x(1))*(t-x(2))*(t-x(4))*(t-x(5)))/((x(3)-x(1))*(x(3)-x(2))*(x(3)-x(4))*(x(3)-x(5)));
l4 = @(t) ((t-x(1))*(t-x(2))*(t-x(3))*(t-x(5)))/((x(4)-x(1))*(x(4)-x(2))*(x(4)-x(3))*(x(4)-x(5)));
l5 = @(t) ((t-x(1))*(t-x(2))*(t-x(3))*(t-x(4)))/((x(5)-x(1))*(x(5)-x(2))*(x(5)-x(3))*(x(5)-x(4)));

% The interpolating polynomial uses the interpolated function values as the
% coefficients. One nice feature!
l = @(t) y(1)*l1(t)+y(2)*l2(t)+y(3)*l3(t)+y(4)*l4(t)+y(5)*l5(t);

% Evaluate at a new point and compare to Lagrange
l(90)

% The Newton polynomials are a sequence of successively higher-order
% polynomials that are defined recursively and are products of factors
% centered at the interpolation points. To solve for the coefficients, a
% lower-triangular system is solved

p1 = @(t) 1;
p2 = @(t) p1(t)*(t-x(1));
p3 = @(t) p2(t)*(t-x(2));
p4 = @(t) p3(t)*(t-x(3));
p5 = @(t) p4(t)*(t-x(4));

% Apply the above functions to the interpolation points and solve the
% system for the coefficients
A = [arrayfun(p1, x); arrayfun(p2, x); arrayfun(p3, x); arrayfun(p4, x); arrayfun(p5, x)]';
c = A\y';
p = @(t) c(1)*p1(t)+c(2)*p2(t)+c(3)*p3(t)+c(4)*p4(t)+c(5)*p5(t);

% Evaluate at a new point and compare to Lagrange
p(90)


