% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 4/8/22
%
% This script graphs the error of the Chebyshev and Legendre least-squares
% approximation to sin(pi/2*x) on the interval [-1,1]

syms x
f = @(x) sin(pi*x/2);
a = -1;
b = 1;
order = 4;

t = linspace(-1, 1);
cheb = least_squares_cheb(f, x, a, b);
leg = least_squares_leg(f, x, a, b);
error_cheb = abs(f(t)-cheb(t));
error_leg = abs(f(t)-leg(t));

plot(t, error_cheb)
hold on
plot(t, error_leg)
grid on
legend('Chebyshev', 'Legendre')
xlabel('Interval [-1,1]');
ylabel('Error');
title('Chebyshev and Legendre Least Squares Approximation');

function [p] = least_squares_cheb(f, x, a, b)
    cheb1 = 1;
    cheb2 = x;
    cheb3 = 2*x^2-1;
    cheb4 = 4*x^3-3*x;
    chebs = [cheb1, cheb2, cheb3, cheb4];
    num_chebs = 4;
    weight = @(x) (1-x^2)^(-1/2);
   
    coeffs = sym(zeros(1, num_chebs));
    p = sym(0);
    for i=1:num_chebs
        coeffs(i) = int(sym(f)*sym(weight)*chebs(i), x, a, b)/int(sym(weight)*chebs(i)*chebs(i), x, a, b);
        p = p + coeffs(i)*chebs(i);
    end
    
    p = matlabFunction(p);
end

function [p] = least_squares_leg(f, x, a, b)
    leg_1 = 1;
    leg_2 = x;
    leg_3 = (3*x^2-1)/2;
    leg_4 = (5*x^3-3*x)/2;
    legs = [leg_1, leg_2, leg_3, leg_4];
    weight = @(x) 1;
    num_legs = 4;
    
    coeffs = sym(zeros(1, num_legs));
    p = sym(0);
    for i=1:num_legs
        coeffs(i) = int(sym(f)*sym(weight)*legs(i), x, a, b)/int(sym(weight)*legs(i)*legs(i), x, a, b);
        p = p + coeffs(i)*legs(i);
    end
    
    p = matlabFunction(p);
end

