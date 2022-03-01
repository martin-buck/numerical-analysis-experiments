% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 2/21/22
%
% This script solves a non-linear system using Newton's Method and verifies
% that the convergence is quadratic by plotting the ratio of errors

% It can be verified analytically that the solution this will converge to
% is [1/sqrt(2), 1/sqrt(2)]
x = [1/sqrt(2); 1/sqrt(2)];
xk = [0; 1];
% This array will store the ratio of errors to verify quadratic convergence
ek_arr = [];
ek = norm(xk-x, 2);

disp(xk)
f1 = @(x, y) x^2 - y^2;
f2 = @(x, y) 2*x*y - 1;

Jf11 = @(x) 2*x;
Jf12 = @(y) -2*y;
Jf21 = @(y) 2*y;
Jf22 = @(x) 2*x;

Jxk = [Jf11(xk(1)), Jf12(xk(2)); Jf21(xk(2)), Jf22(xk(1))];
f = [f1(xk(1), xk(2)); f2(xk(1), xk(2))];
xk1 = Jxk\(-f) + xk;
ek1 = norm(xk1-x, 2);
count = 1;
format long
disp(xk1)

while(ek1 > 10^-8)
    xk = xk1;
    ek = norm(xk-x, 2);
    Jxk = [Jf11(xk(1)), Jf12(xk(2)); Jf21(xk(2)), Jf22(xk(1))];
    f = [f1(xk(1), xk(2)); f2(xk(1), xk(2))];
    xk1 = Jxk\(-f) + xk;
    ek1 = norm(xk1-x, 2);
    ek_arr(count) = ek1/ek^2;
    count = count + 1;
    disp(xk1)
end

scatter(1:count-1, ek_arr, 'filled')
grid on
ylabel('Error Ratio $\frac{||e_k||}{||e_{k-1}||^2}$', 'Interpreter', 'latex')
xlabel('Iteration Number', 'Interpreter', 'latex')
title(strcat('Error Ratio $\frac{||e_k||}{||e_{k-1}||^2}$ vs. Iteration Number'), 'Interpreter', 'latex')
hold on
dx = .01; dy = .01;
text(1:count-1+dx, ek_arr+dy, cellstr(num2str(ek_arr')))