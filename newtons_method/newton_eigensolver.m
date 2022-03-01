% Martin Buck
% Tufts University - Math 225 - Numerical Analysis
% 2/21/22
%
% This script computes an eigenpair iteratively using Newton's Method

A = 10*rand(100);
x0 = 10*rand(100,1);
x0 = x0/norm(x0, 2);

E = eig(A);
[x, l, count] = newton_eigenpair(A, x0);

function [x, l, count] = newton_eigenpair(A, x0)
    [m, n] = size(A);
    
    lk = x0'*A*x0;
    xk = x0;
    B = [A-lk*eye(m,n), -xk; 2*xk', 0];
    b = -[(A*xk)-(lk*xk); (xk'*xk)-1];
    x_temp = B\b;
    lk1 = lk + x_temp(end);
    xk1 = xk + x_temp(1:end-1);
    count = 0;
    
    while(max(abs(xk1-xk)) > (10^-8))
        xk = xk1;
        lk = lk1;
        B = [A-lk*eye(m,n), -xk; 2*xk', 0];
        b = -[(A*xk)-(lk*xk); (xk'*xk)-1];
        x_temp = B\b;
        lk1 = lk + x_temp(end);
        xk1 = xk + x_temp(1:end-1);
        count = count + 1;
    end
    
    x = xk1;
    l = lk1;
end


    