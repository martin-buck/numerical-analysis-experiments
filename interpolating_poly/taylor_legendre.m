% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 4/2/22
%
% This script plots the error for a Taylor polynomial expansion of ln(x) on
% the interval [1,2]. It also plots the error of the Legendre polynomial
% approximation on the same interval. This can be computed analytically
% using an orthogonal projection

% First-order Taylor polynomial
T = @(x) log(3/2)+(2/3)*(x-3/2);

% Least squares approximation by a first-order polynomial using the
% Legendre basis
L = @(x) (log(4)-1) - 3*(log(16)-3)*(x-3/2);

t = linspace(1,2,50);

plot(t, abs(T(t)-log(t)));
hold on
plot(t, abs(L(t)-log(t)));
legend('First Order Taylor', 'First Order Legendre')
grid on
xlabel('x')
ylabel('Approximation Error')
title('Taylor and Legendre Polynomial Approximation to Logarithm')